#!/usr/bin/env julia

using Base.Meta

iscall(e::Expr, name::Symbol) = isexpr(e,:call) && e.args[1] == name
iscall(e::Expr, names::Vector) = isexpr(e,:call) && e.args[1] in names
iscall(x, y) = false

map_expr(f::Function, e::Expr) = f(Expr(e.head, [map_expr(f,a) for a in e.args]...))
map_expr(f::Function, x) = x

map_nonexpr(f::Function, e::Expr) = Expr(e.head, [map_nonexpr(f,a) for a in e.args]...)
map_nonexpr(f::Function, x) = f(x)

function mapreduce_expr(f::Function, op::Function, e::Expr)
    s = f(e)
    for a in e.args
        if isa(a,Expr)
            s = op(s, mapreduce_expr(f, op, a))
        end
    end
    s
end

filter_expr(f::Function, e::Expr) = Expr(e.head, [filter_expr(f,x) for x in filter(f,e.args)]...)
filter_expr(f::Function, x) = x

function unravel_data(e::Expr)
    if e.head == :macrocall && e.args[1] == Symbol("@DATA")
        vars = []
        vals = []
        @assert e.args[2].head == :tuple
        b = true
        for x in e.args[2].args
            if isa(x,Expr) && x.head == :call
                @assert x.args[1] == :/
                push!(vars, x.args[2])
                push!(vals, x.args[3])
                b = !b
            else
                b ? push!(vars,x) : push!(vals,x)
            end
        end
        return Expr(:block, map((a,b)->Expr(:(=),a,b), vars, vals)...)
    end
    e
end

mutable struct FortranRoutine
    filename::String
    ast::Expr
    name::Symbol
    args::Vector{Symbol}
    isdimargs::Vector{Bool}         # arg is array
    isoutargs::Vector{Bool}         # arg is scalar & output
    issubroutine::Bool              # function or subroutine?
    calls::Set{Symbol}
    types::Dict{Symbol,Symbol}
    dims::Dict{Symbol,Vector{Any}}
    locals::Set{Symbol}
    globals::Vector{Expr}
end

 GLOBALS = Set([:D1MACH,:I1MACH])
 INTRINSICS = Set([:COMPLEX,:DABS,:DATAN,:DBLE,:DCOS,:DCOSH,:DEXP,
    :DLOG,:DMAX1,:DMIN1,:DSIGN,:DSIN,:DSINH,:DSQRT,:FLOAT,:IABS,:INT,
    :MAX0,:MIN0,:MOD,:SNGL])

function main()

    routines = Dict{Symbol,FortranRoutine}()

    for filename in filter(s->endswith(s,".jl0"), readdir())
        src = read(filename, String)
        ast = Base.Meta.parse(src)

        @assert ast.head == :function
        name = ast.args[1].args[1]
        args = ast.args[1].args[2:end]

        # strip line numbers
        ast = filter_expr(e->!(isa(e,LineNumberNode) || isexpr(e,:line)), ast)

        # D1MACH(i) --> D1MACHi
        ast = map_expr(
                e->(!(iscall(e,:D1MACH) || iscall(e,:I1MACH)) ? e : Symbol(string(e.args[1], e.args[2]::Int))),
                ast)

        # X(i) = block Y end --> X[i] = Y
        ast = map_expr(
                e->(!(isexpr(e,:(=)) &&
                      isexpr(e.args[1],:call) &&
                      isexpr(e.args[2],:block) &&
                      length(e.args[2].args) == 1) ? e :
                    Expr(:(=), Expr(:ref, e.args[1].args...), e.args[2].args[1])),
                ast)

        # extract type & dimension info

        function extract_macro(name)
            # get macro args
            comma_args(e) = isexpr(e,:tuple) ? e.args : [e]
            margs = mapreduce_expr(
                e->(isexpr(e,:macrocall) && e.args[1] == name ?
                    Set(comma_args(e.args[2])) : Set()),
                union!, ast)
            # delete macro expr
            ast = filter_expr(e->!(isexpr(e,:macrocall) && e.args[1] == name), ast)
            margs
        end

        vartypes = Dict{Symbol,Symbol}()
        for (s,t) in ((Symbol("@DOUBLE_PRECISION"), :Float64),
                      (Symbol("@INTEGER"), :Int32))
            merge!(vartypes, Dict{Symbol,Symbol}(v => t for v in extract_macro(s)))
        end
        vardims = Dict(e.args[1] => e.args[2:end] for e in extract_macro(Symbol("@DIMENSION")))

        # add types to signature
        ts = [get(vartypes,v,:Int32) for v in args]
        ast.args[1].args[2:end] = [ haskey(vardims,v) ? :($v::AbstractArray{$t}) : :($v::$t) for (v,t) in zip(args,ts) ]

        # X(i) --> X[i] when X is a dimensioned variable
        dv = collect(keys(vardims))
        ast = map_expr(
                e->(!iscall(e,dv) ? e : Expr(:ref, e.args...)),
                ast)

        # get calls
        calls = mapreduce_expr(
            e->(isexpr(e,:call) ?  Set([e.args[1]]) : Set(Symbol[])),
            union!, ast)

        # fix sloppy ZABS signature: ZABS(ZR,ZI) --> ZABS(Z)
        if name == :ZABS
            args = [:Z]
            ast.args[1].args = [:ZABS, :(Z::Complex128)]
            pushfirst!(ast.args[2].args, :(ZI = imag(Z)))
            pushfirst!(ast.args[2].args, :(ZR = real(Z)))
        end

        # @DATA x, y / 1, 2 --> x = 1; y = 2
        ast = map_expr(unravel_data, ast)

        # x/y --> div(x,y) if x,y are ints
        ast = map_expr(
            function (e)
                if iscall(e,:/) &&
                   (get(vartypes,e.args[2],nothing) == :Int32 || isa(e.args[2],Int)) &&
                   (get(vartypes,e.args[3],nothing) == :Int32 || isa(e.args[3],Int))
                    @assert length(e.args) == 3
                    return Expr(:call, :div, e.args[2], e.args[3])
                end
                e
            end, ast)

        # convert all int literals to int32
        ast = map_nonexpr(x->(isa(x,Int) ? :(int32($x)) : x), ast)

        # get scalar assignments
        modifies = mapreduce_expr(
            e->(isexpr(e,:(=)) && isa(e.args[1],Symbol) ?
                Set([e.args[1]]) : Set(Symbol[])),
            union!, ast)

        isdimargs = [haskey(vardims,a) for a in args]
        isoutargs = [a in modifies for a in args]
        issubroutine = !(name in modifies)

        routines[name] = FortranRoutine(filename, ast, name, args, isdimargs, isoutargs,
                                        issubroutine, calls, vartypes, vardims, Set{Symbol}(),
                                        Expr[])
    end

    # figure out and initialize local vars

    for (name,routine) in routines
        for v in keys(routine.types)
            if !(v in GLOBALS || v in INTRINSICS ||
                 haskey(routines,v) || v in routine.args)
                push!(routine.locals, v)
            end
        end

        for v in sort(collect(routine.locals), rev=true)
            t = routine.types[v]
            if haskey(routine.dims,v)
                d = routine.dims[v]
                gv = Symbol("_$(name)_$v")
                push!(routine.globals, :( $gv = Array($t,$(d...))))
                e = :( $v = $gv)
            else
                e = :($v::$t = zero($t))
            end
            pushfirst!(routine.ast.args[2].args, e)
        end

        # function return value
        if !routine.issubroutine
            v = routine.name
            t = routine.types[v]
            pushfirst!(routine.ast.args[2].args, :($v::$t = zero($t)))
        end
    end

    # propagate isoutargs across calls
    # (turns out this wasn't necessary)

    for _ = 1:0 # blindly iterate a few times to converge
        for (name,routine) in routines
            for i in find(!routine.isoutargs)
                arg = routine.args[i]
                b = mapreduce_expr(
                        function (e)
                            if isexpr(e,:call) && haskey(routines,e.args[1])
                                r = routines[e.args[1]]
                                idx = findin([arg], e.args[2:end])
                                return any(r.isoutargs[idx])
                            end
                            false
                        end, |, routine.ast)
                b && println("$name $arg is modified")
                routine.isoutargs[i] = b
            end
        end
    end

    # pull outargs out of subroutine calls

    for (name,routine) in routines
        any(routine.isoutargs) || continue
        routine.issubroutine || continue

        # call F(A,B) --> B = F(A,B) if B is an outarg
        for (name2,routine2) in routines
            (name in routine2.calls) || continue
            ioa = routine.isoutargs
            routine2.ast.args[2] = map_expr(
                e->(!iscall(e,name) ? e :
                    Expr(:(=), Expr(:tuple, e.args[2:end][ioa]...), e)),
                routine2.ast.args[2])

            # hack to replace C[I,J] --> slice(X,I:end,J) in subroutine calls
            routine2.ast = map_expr(
                function (e)
                    iscall(e,name) || return e
                    newargs = []
                    for (i,a) in enumerate(e.args[2:end])
                        if isexpr(a,:ref) && routine.isdimargs[i]
                            @assert length(a.args) == 3
                            a = Expr(:call, :slice,
                                     a.args[1],
                                     Expr(:(:), a.args[2], routine2.dims[a.args[1]][1]),
                                     Expr(:call, :int, a.args[3]))
                        end
                        push!(newargs, a)
                    end
                    Expr(:call, name, newargs...)
                end, routine2.ast)
        end
    end

    # fix return statements

    for (name,routine) in routines
        outargs = routine.args[routine.isoutargs]
        if !routine.issubroutine
            #pushfirst!(outargs, name)
            outargs = [name]
        end
        length(outargs) > 0 ||continue

        # add outargs to returns
        outexpr = Expr(:return, length(outargs) == 1 ? outargs[1] : Expr(:tuple, outargs...))
        routine.ast = map_expr(e->(isexpr(e,:return) ? outexpr : e), routine.ast)

        # mangle function return variable
        if !routine.issubroutine
            namevar = Symbol("__$(name)__")
            routine.ast.args[2] = map_nonexpr(
                s->(s != name ? s : namevar),
                routine.ast.args[2])
        end
    end

    # write .jl files

    for (_,routine) in routines
        fn = splitext(routine.filename)[1] * ".jl"
        open(fn,"w") do f
            for g in routine.globals
                println(f, g)
            end
            println(f, routine.ast)
        end
    end
end

if !isinteractive()
    main()
end

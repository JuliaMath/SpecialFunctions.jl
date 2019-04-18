# Building OpenSpecFun from scratch

using BinDeps
using BinDeps: libdir, srcdir, includedir, depsdir, builddir
using Base.Math: libm
import BinaryProvider
using Libdl

BinDeps.@setup
const OSF_VERS = v"0.5.3"
openspecfun = library_dependency("libopenspecfun")

# If Julia is built with OpenLibm, we want to build OpenSpecFun with it as well.
# Unfortunately this requires a fair bit more work, as we need to link to the .so
# and to include the headers, which aren't readily available.
if libm == "libopenlibm"
    const OLM_VERS = v"0.5.4"
    const OLM_HASH = "9a8ae1d17825a4a6a4c013d36a7f4348b27c47eedb6549c521ecc9c79d021c13"
    const OLM_URL = "https://github.com/JuliaLang/openlibm/archive/v$OLM_VERS.tar.gz"
    use_openlibm = true

    if !isdir(libdir(openspecfun))
        mkpath(libdir(openspecfun))
    end

    openlibm_so = Libdl.dlpath(libm)

    # Copy over the OpenLibm .so
    for lib in readdir(dirname(openlibm_so))
        startswith(lib, "libopenlibm") || continue
        cp(joinpath(dirname(openlibm_so), lib), joinpath(libdir(openspecfun), lib),
           force=true, follow_symlinks=false)
    end

    # Grab and unpack the tarball so we can get the header filesa
    BinaryProvider.download_verify_unpack(OLM_URL, OLM_HASH, srcdir(openspecfun); verbose = true)

    # Copy over all of the OpenLibm headers
    cp(joinpath(srcdir(openspecfun), "openlibm-$OLM_VERS"), joinpath(srcdir(openspecfun), "openlibm"); force=true)
    openlibm_src = joinpath(srcdir(openspecfun), "openlibm")
    openlibm_include = joinpath(includedir(openspecfun), "openlibm")
    if !isdir(openlibm_include)
        mkpath(openlibm_include)
    end
    for f in readdir(joinpath(openlibm_src, "include"))
        cp(joinpath(openlibm_src, "include", f), joinpath(openlibm_include, f),
           force=true, follow_symlinks=true)
    end
    for f in readdir(joinpath(openlibm_src, "src"))
        if endswith(f, ".h")
            cp(joinpath(openlibm_src, "src", f), joinpath(openlibm_include, f),
               force=true, follow_symlinks=true)
        end
    end
else
    use_openlibm = false
end

fc = "gfortran"

# macOS has precompiled binaries, so it's just FreeBSD that should default to Clang
if Sys.KERNEL === :FreeBSD || Sys.KERNEL == :Darwin
    cc = "clang"
    use_clang = true
else
    cc = "gcc"
    use_clang = false
end

if Sys.ARCH in [:i386, :i387, :i486, :i586, :i686]
    cc *= " -m32"
    fc *= " -m32"
elseif Sys.ARCH === :x86_64
    cc *= " -m64"
    fc *= " -m64"
end

extra_ld = Sys.isapple() ? "" : "-Wl,-rpath,'\$\$ORIGIN' -Wl,-z,origin"

provides(Sources, URI("https://github.com/JuliaLang/openspecfun/archive/v$OSF_VERS.tar.gz"),
         openspecfun, unpacked_dir="openspecfun-$OSF_VERS")

provides(BuildProcess,
    (@build_steps begin
        GetSources(openspecfun)
        CreateDirectory(builddir(openspecfun))
        @build_steps begin
            ChangeDirectory(joinpath(srcdir(openspecfun), "openspecfun-$OSF_VERS"))
            FileRule(joinpath(libdir(openspecfun), "libopenspecfun." * Libdl.dlext),
                @build_steps begin
                    CreateDirectory(libdir(openspecfun))
                    ```
                    $MAKE_CMD install
                        ARCH="$(Sys.ARCH)"
                        CC="$cc"
                        FC="$fc"
                        USECLANG=$(Int(use_clang))
                        USEGCC=$(Int(!use_clang))
                        USE_OPENLIBM=$(Int(use_openlibm))
                        CFLAGS="-O3 -std=c99"
                        FFLAGS="-O2 -fPIC"
                        LDFLAGS="-L$(libdir(openspecfun)) $extra_ld"
                        DESTDIR=""
                        prefix=""
                        libdir="$(libdir(openspecfun))"
                        includedir="$(includedir(openspecfun))"
                        O=
                    ```
                end)
        end
    end), openspecfun)

BinDeps.@install Dict(:libopenspecfun => :openspecfun)

# BinDeps doesn't give us a `check_deps` function in the generated deps.jl file like
# BinaryProvider does. Instead, it uses a macro and checks immediately whether the
# library can be loaded. So let's just fake one.
depsjl = joinpath(depsdir(openspecfun), "deps.jl")
@assert isfile(depsjl)
if !any(line->occursin("check_deps", line), eachline(depsjl))
    open(depsjl, "a") do fh
        println(fh, """
            # NOTE: This function is a compatibility shim for BinaryProvider-based builds
            check_deps() = nothing
            """)
    end
end

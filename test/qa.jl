@testset "Aqua" begin
    Aqua.test_all(SpecialFunctions)
end

@testset "ExplicitImports" begin
    # No implicit imports (`using XY`)
    @test ExplicitImports.check_no_implicit_imports(SpecialFunctions) === nothing

    # All explicit imports (`using XY: Z`) are loaded via their owners
    @test ExplicitImports.check_all_explicit_imports_via_owners(
        SpecialFunctions;
        ignore = (
            # Ref https://github.com/JuliaTesting/ExplicitImports.jl/issues/92
            :invπ, # SpecialFunctions
            :sqrtπ, # SpecialFunctions
        ),
    ) === nothing

    # Limit explicit imports (`using XY: Z`) of non-public names to a minimum
    @test ExplicitImports.check_all_explicit_imports_are_public(
        SpecialFunctions;
        ignore = (
            :MPFRRoundingMode, # Base.MPFR
            :ROUNDING_MODE, # Base.MPFR
            :nan_dom_err, # Base.Math
            # Ref https://github.com/JuliaTesting/ExplicitImports.jl/issues/92
            :invπ, # SpecialFunctions
            :sqrtπ, # SpecialFunctions
        ),
    ) === nothing

    # No explicit imports (`using XY: Z`) that are not used
    @test ExplicitImports.check_no_stale_explicit_imports(SpecialFunctions) === nothing

    # Nothing is accessed via modules other than its owner
    @test ExplicitImports.check_all_qualified_accesses_via_owners(SpecialFunctions) === nothing

    # Limit accesses of non-public names to a minimum
    @test ExplicitImports.check_all_qualified_accesses_are_public(
        SpecialFunctions;
        ignore = (
            :IEEEFloat, # Base
            :MPFR, # Base
            :MPFRRoundingMode, # Base.MPFR
            :ROUNDING_MODE, # Base.MPFR
            :_fact_table64, # Base
            :version, # Base.MPFR
            Symbol("@nif"), # Base
            (VERSION < v"1.11" ? (:depwarn,) : ())..., # Base
        ),
    ) === nothing

    # No self-qualified accesses
    @test ExplicitImports.check_no_self_qualified_accesses(SpecialFunctions) === nothing
end

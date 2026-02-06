using JET: JET

@testset "package" begin
    # Check that there are no undefined global references and undefined field accesses
    JET.test_package(SpecialFunctions; target_defined_modules = true, mode = :typo)

    # Analyze methods based on their declared signature
    JET.report_package(SpecialFunctions; target_defined_modules = true)
end

@testset "logabsgamma" begin
    # issue #502
    JET.@test_call logabsgamma(1.0)
    JET.@test_opt logabsgamma(1.0)
    JET.@test_call logabsgamma(1f0)
    JET.@test_opt logabsgamma(1f0)
end

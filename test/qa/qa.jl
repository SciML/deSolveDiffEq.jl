using deSolveDiffEq
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(deSolveDiffEq)
end

@testset "JET" begin
    JET.test_package(deSolveDiffEq; target_defined_modules = true)
end

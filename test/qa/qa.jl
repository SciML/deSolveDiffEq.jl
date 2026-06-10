using deSolveDiffEq
using Aqua
using JET
using Test

@testset "Aqua" begin
    # persistent_tasks disabled: RCall's __init__ calls `rimport`, which `eval`s
    # into a module and breaks incremental precompilation of the child package
    # Aqua spawns, so the check cannot pass until RCall init is precompile-safe.
    Aqua.test_all(deSolveDiffEq; persistent_tasks = false)
    @test_broken false  # Aqua persistent_tasks: RCall __init__ eval breaks incremental precompile — see https://github.com/SciML/deSolveDiffEq.jl/issues/51
end

@testset "JET" begin
    JET.test_package(deSolveDiffEq; target_defined_modules = true)
end

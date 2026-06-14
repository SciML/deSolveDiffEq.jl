using deSolveDiffEq
using AllocCheck
using Test

# Allocation tests for pure Julia helper functions
# Note: The main solve function allocates due to R interop (RCall),
# so we only test the Julia-side helper functions
@testset "AllocCheck - Algorithm Name Functions" begin
    # Test that algname functions are allocation-free (compile-time constants)
    algs = (
        deSolveDiffEq.lsoda(), deSolveDiffEq.lsode(), deSolveDiffEq.lsodes(),
        deSolveDiffEq.lsodar(), deSolveDiffEq.vode(), deSolveDiffEq.daspk(),
        deSolveDiffEq.euler(), deSolveDiffEq.rk4(), deSolveDiffEq.ode23(),
        deSolveDiffEq.ode45(), deSolveDiffEq.radau(), deSolveDiffEq.bdf(),
        deSolveDiffEq.bdf_d(), deSolveDiffEq.adams(), deSolveDiffEq.impAdams(),
        deSolveDiffEq.impAdams_d(), deSolveDiffEq.iteration(),
    )

    for alg in algs
        # Warmup
        deSolveDiffEq.algname(alg)
        # Check that algname doesn't allocate
        allocs = @allocated deSolveDiffEq.algname(alg)
        @test allocs == 0
    end
end

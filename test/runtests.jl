using deSolveDiffEq
using Test

function lorenz(u, p, t)
    du1 = 10.0(u[2] - u[1])
    du2 = u[1] * (28.0 - u[3]) - u[2]
    du3 = u[1] * u[2] - (8 / 3) * u[3]
    return [du1, du2, du3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, deSolveDiffEq.lsoda())
sol = solve(prob, deSolveDiffEq.lsode())
sol = solve(prob, deSolveDiffEq.lsodes())
sol = solve(prob, deSolveDiffEq.lsodar())
sol = solve(prob, deSolveDiffEq.vode())
sol = solve(prob, deSolveDiffEq.daspk())
sol = solve(prob, deSolveDiffEq.euler())
sol = solve(prob, deSolveDiffEq.rk4())
sol = solve(prob, deSolveDiffEq.ode23())
sol = solve(prob, deSolveDiffEq.ode45())
sol = solve(prob, deSolveDiffEq.radau())
sol = solve(prob, deSolveDiffEq.bdf())
sol = solve(prob, deSolveDiffEq.bdf_d())
sol = solve(prob, deSolveDiffEq.adams())
sol = solve(prob, deSolveDiffEq.impAdams())
sol = solve(prob, deSolveDiffEq.impAdams_d())
#sol = solve(prob,deSolveDiffEq.iteration())

sol = solve(prob, deSolveDiffEq.lsoda(), saveat = 0.1)
#using Plots; plot(sol,vars=(1,2,3))

# Allocation tests for pure Julia helper functions
# Note: The main solve function allocates due to R interop (RCall),
# so we only test the Julia-side helper functions
if get(ENV, "GROUP", "all") == "all" || get(ENV, "GROUP", "all") == "nopre"
    using AllocCheck

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
end

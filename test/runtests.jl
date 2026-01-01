using deSolveDiffEq
using Test

@testset "deSolveDiffEq.jl" begin
    @testset "Algorithm types" begin
        # Test type hierarchy
        @test lsoda() isa deSolveAlgorithm
        @test lsoda() isa DiffEqBase.AbstractODEAlgorithm
        @test lsode() isa deSolveAlgorithm
        @test lsodes() isa deSolveAlgorithm
        @test lsodar() isa deSolveAlgorithm
        @test vode() isa deSolveAlgorithm
        @test daspk() isa deSolveAlgorithm
        @test euler() isa deSolveAlgorithm
        @test rk4() isa deSolveAlgorithm
        @test ode23() isa deSolveAlgorithm
        @test ode45() isa deSolveAlgorithm
        @test radau() isa deSolveAlgorithm
        @test bdf() isa deSolveAlgorithm
        @test bdf_d() isa deSolveAlgorithm
        @test adams() isa deSolveAlgorithm
        @test impAdams() isa deSolveAlgorithm
        @test impAdams_d() isa deSolveAlgorithm
        @test iteration() isa deSolveAlgorithm
    end

    @testset "ODE solving" begin
        function lorenz(u, p, t)
            du1 = 10.0(u[2]-u[1])
            du2 = u[1]*(28.0-u[3]) - u[2]
            du3 = u[1]*u[2] - (8/3)*u[3]
            [du1, du2, du3]
        end
        u0 = [1.0; 0.0; 0.0]
        tspan = (0.0, 100.0)
        prob = ODEProblem(lorenz, u0, tspan)

        sol = solve(prob, lsoda())
        @test sol.retcode == ReturnCode.Default || sol.retcode == ReturnCode.Success
        sol = solve(prob, lsode())
        sol = solve(prob, lsodes())
        sol = solve(prob, lsodar())
        sol = solve(prob, vode())
        sol = solve(prob, daspk())
        sol = solve(prob, euler())
        sol = solve(prob, rk4())
        sol = solve(prob, ode23())
        sol = solve(prob, ode45())
        sol = solve(prob, radau())
        sol = solve(prob, bdf())
        sol = solve(prob, bdf_d())
        sol = solve(prob, adams())
        sol = solve(prob, impAdams())
        sol = solve(prob, impAdams_d())
        #sol = solve(prob, iteration())

        sol = solve(prob, lsoda(), saveat = 0.1)
        @test length(sol.t) > 2
    end
end

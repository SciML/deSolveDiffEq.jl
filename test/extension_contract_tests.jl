using deSolveDiffEq
using SciMLBase
using Test

module ExternalSciMLBaseExtension

    using SciMLBase

    struct MockAlgorithm <: SciMLBase.AbstractODEAlgorithm end

    SciMLBase.__solve(
        ::SciMLBase.AbstractODEProblem,
        ::MockAlgorithm,
        args...;
        kwargs...,
    ) = :external_extension

end

@testset "Public solver extension contract" begin
    f(u, p, t) = u
    prob = SciMLBase.ODEProblem(f, 1.0, (0.0, 1.0))

    @test SciMLBase.solve(prob, ExternalSciMLBaseExtension.MockAlgorithm()) ===
        :external_extension
end

@testset "Public algorithms" begin
    algorithms = (
        lsoda, lsode, lsodes, lsodar, vode, daspk, euler, rk4, ode23, ode45, radau,
        bdf, bdf_d, adams, impAdams, impAdams_d, iteration,
    )

    @test all(algorithm() isa deSolveAlgorithm for algorithm in algorithms)
end

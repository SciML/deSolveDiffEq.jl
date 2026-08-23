using SciMLTesting, deSolveDiffEq, Test

# ExplicitImports can only check an extension module that actually exists, and an
# extension module only exists once its triggers are loaded. Without this `using`
# the ext/ sources are never scanned by QA at all. RCall needs a working R
# installation, which CI provides through the `r-base-dev` apt package.
using RCall

# ExplicitImports silently skips an extension that fails to load, so assert the
# extension modules actually exist rather than trusting a green run_qa.
@testset "Extensions loaded" begin
    for ext in (:deSolveDiffEqRCallExt,)
        @test Base.get_extension(deSolveDiffEq, ext) !== nothing
    end
end

# The SciML common interface deSolveDiffEq deliberately reexports so that
# `using deSolveDiffEq` is enough to build an ODE problem, solve it, and inspect the
# result. Owned and documented upstream; kept in sync with the reexport `export` block
# in src/deSolveDiffEq.jl.
const REEXPORTS = (
    :DEStats, :EnsembleAnalysis, :EnsembleDistributed, :EnsembleProblem, :EnsembleSerial,
    :EnsembleSolution, :EnsembleSplitThreads, :EnsembleSummary, :EnsembleThreads,
    :NullParameters, :ODEFunction, :ODEProblem, :ODESolution, :ReturnCode, :remake,
    :solve, :successful_retcode,
)

run_qa(
    deSolveDiffEq;
    reexports_allow = REEXPORTS,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                # deSolveDiffEq's own internal, reached from its own extension.
                # ExplicitImports treats an extension as a separate module, so this
                # reads as a non-public cross-module access even though it never
                # leaves the package.
                :r_adapter,
            ),
        ),
    ),
)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using deSolveDiffEq`, so
    # the allow-list cannot drift into approving names the package no longer provides.
    # `isdefined(@__MODULE__, ...)` tests the property directly: this file's
    # `using deSolveDiffEq` is what has to bring the name into scope.
    @testset "$name" for name in REEXPORTS
        @test name in names(deSolveDiffEq)
        @test isdefined(@__MODULE__, name)
    end
end

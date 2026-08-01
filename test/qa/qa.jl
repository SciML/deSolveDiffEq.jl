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

run_qa(
    deSolveDiffEq;
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

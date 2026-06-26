using SciMLTesting, deSolveDiffEq, Test
using JET

run_qa(
    deSolveDiffEq;
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    # persistent_tasks: RCall's __init__ calls `rimport`, which `eval`s into a
    # module and breaks incremental precompilation of the child package Aqua
    # spawns, so the check cannot pass until RCall init is precompile-safe.
    # https://github.com/SciML/deSolveDiffEq.jl/issues/51
    aqua_broken = (:persistent_tasks,),
    ei_kwargs = (;
        # AbstractODEAlgorithm/AbstractODEProblem/__solve/build_solution are owned
        # by SciMLBase but reached through DiffEqBase (this package's declared dep,
        # not SciMLBase); they are the canonical DiffEqBase solve interface this
        # package subtypes/extends/calls.
        all_qualified_accesses_via_owners = (;
            ignore = (:AbstractODEAlgorithm, :AbstractODEProblem, :__solve, :build_solution),
        ),
        # Same four names accessed as `DiffEqBase.X` are non-public in DiffEqBase
        # (owned by SciMLBase, re-exported non-public through DiffEqBase).
        all_qualified_accesses_are_public = (;
            ignore = (:AbstractODEAlgorithm, :AbstractODEProblem, :__solve, :build_solution),
        ),
    ),
)

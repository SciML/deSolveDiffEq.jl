using SciMLTesting, deSolveDiffEq, Test
using JET

run_qa(
    deSolveDiffEq;
    api_docs_kwargs = (;
        rendered = true,
        rendered_ignore = Tuple(names(deSolveDiffEq.DiffEqBase)),
    ),
    explicit_imports = true,
    jet_kwargs = (; target_defined_modules = true),
    # persistent_tasks: RCall's __init__ calls `rimport`, which `eval`s into a
    # module and breaks incremental precompilation of the child package Aqua
    # spawns, so the check cannot pass until RCall init is precompile-safe.
    # https://github.com/SciML/deSolveDiffEq.jl/issues/51
    aqua_broken = (:persistent_tasks,),
    ei_kwargs = (;
        # `__solve` is the SciMLBase-owned solve-interface hook this package
        # extends as `SciMLBase.__solve`. It remains non-public in SciMLBase (and
        # in DiffEqBase), so the public-API check flags the qualified access even
        # though it is reached through its owner.
        all_qualified_accesses_are_public = (;
            ignore = (:__solve,),
        ),
    ),
)

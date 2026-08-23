# API

The deSolveDiffEq-owned algorithms are documented on the [home page](index.md).

## Reexported SciML common interface

`using deSolveDiffEq` also brings in the parts of the SciML common interface needed to
build an ODE problem, solve it, and inspect the result, so they do not have to be
imported separately. deSolveDiffEq does not define these names -- they are owned and
documented by [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/), and that is where
their documentation lives:

  - Problems: `ODEProblem`, `EnsembleProblem`
  - Functions: `ODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`, `DEStats`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - `NullParameters`

`RCall` still has to be loaded separately -- it is a weak dependency, and the
`deSolveDiffEqRCallExt` extension is what wires up the R-backed solvers.

Anything else from SciMLBase must be imported from SciMLBase directly. Three groups are
deliberately absent:

  - **DAE, SDE, DDE and every other non-ODE problem type.** `SciMLBase.__solve` is
    defined here only for `AbstractODEProblem`.
  - **Callbacks.** deSolveDiffEq passes no callback through to R, so
    `ContinuousCallback` and friends are not part of its surface.
  - **The integrator interface** (`init`, `step!`, `solve!`, `reinit!`, ...).
    deSolveDiffEq implements `SciMLBase.__solve` only; it has no integrator.

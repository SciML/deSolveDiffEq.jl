# deSolveDiffEq.jl

`deSolveDiffEq.jl` integrates the R [deSolve](https://desolve.r-forge.r-project.org/)
solvers with the SciML solve interface. Load `RCall` before solving so the package's
RCall extension can initialize the R-backed algorithms.

## Algorithms

```@docs
deSolveAlgorithm
lsoda
lsode
lsodes
lsodar
vode
daspk
euler
rk4
ode23
ode45
radau
bdf
bdf_d
adams
impAdams
impAdams_d
iteration
```

```jldoctest; setup = :(using deSolveDiffEq)
julia> deSolveDiffEq.lsoda() isa deSolveDiffEq.deSolveAlgorithm
true
```

## Solver Extension Contract

Solver packages implement the documented `SciMLBase.__solve` interface for their
algorithm types. `deSolveDiffEq` implements that interface for every subtype of
`deSolveAlgorithm`; callers should use `SciMLBase.solve`.

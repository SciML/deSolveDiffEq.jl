using Documenter
using deSolveDiffEq

makedocs(
    sitename = "deSolveDiffEq.jl",
    modules = [deSolveDiffEq],
    checkdocs = :exports,
    doctest = true,
)

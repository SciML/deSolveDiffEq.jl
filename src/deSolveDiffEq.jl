module deSolveDiffEq

using Reexport, RCall
@reexport using DiffEqBase

solver = Ref{Module}()

function __init__()
    solver[] = rimport("deSolve")
end

abstract type deSolveAlgorithm <: DiffEqBase.AbstractODEAlgorithm end

struct lsoda <: deSolveAlgorithm end
struct lsode <: deSolveAlgorithm end
struct lsodes <: deSolveAlgorithm end
struct lsodar <: deSolveAlgorithm end
struct vode <: deSolveAlgorithm end
struct daspk <: deSolveAlgorithm end
struct euler <: deSolveAlgorithm end
struct rk4 <: deSolveAlgorithm end
struct ode23 <: deSolveAlgorithm end
struct ode45 <: deSolveAlgorithm end
struct radau <: deSolveAlgorithm end
struct bdf <: deSolveAlgorithm end
struct bdf_d <: deSolveAlgorithm end
struct adams <: deSolveAlgorithm end
struct impAdams <: deSolveAlgorithm end
struct impAdams_d <: deSolveAlgorithm end
struct iteration <: deSolveAlgorithm end

function DiffEqBase.__solve(
        prob::DiffEqBase.AbstractODEProblem,
        alg::deSolveAlgorithm, timeseries = [], ts = [], ks = [];
        saveat = eltype(prob.tspan)[], timeseries_errors = true,
        reltol = 1e-3, abstol = 1e-6,
        maxiters = 100000,
        kwargs...)
    p = prob.p
    tspan = prob.tspan
    u0 = prob.u0

    if DiffEqBase.isinplace(prob)
        f = function (t, u, __p)
            du = similar(u)
            prob.f(du, u, p, t)
            R"list($du)" # Error message says a list return is required!
        end
    else
        f = function (t, u, __p)
            du = prob.f(u, p, t)
            R"list($du)" # Error message says a list return is required!
        end
    end

    _saveat = isempty(saveat) ? nothing : saveat
    if _saveat isa Array
        __saveat = _saveat
    elseif _saveat isa Number
        __saveat = Array(tspan[1]:_saveat:tspan[2])
    elseif _saveat isa Nothing
        __saveat = [tspan[1], tspan[2]]
    else
        __saveat = Array(_saveat)
    end

    out = rcopy(solver[].ode(times = __saveat, y = u0, func = f,
        method = string(alg)[15:(end - 2)],
        parms = nothing, maxsteps = maxiters,
        rtol = reltol, atol = abstol))

    ts = @view out[:, 1]

    if u0 isa AbstractArray
        timeseries = Vector{typeof(u0)}(undef, length(ts))
        for i in 1:length(ts)
            timeseries[i] = @view out[i, 2:end]
        end
    else
        timeseries = out[i, end]
    end

    DiffEqBase.build_solution(prob, alg, ts, timeseries,
        dense = false,
        timeseries_errors = timeseries_errors)
end

end # module

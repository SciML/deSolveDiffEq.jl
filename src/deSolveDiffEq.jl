module deSolveDiffEq

using Reexport, RCall, PrecompileTools
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

# Compile-time algorithm name extraction to avoid runtime string allocations
algname(::lsoda) = "lsoda"
algname(::lsode) = "lsode"
algname(::lsodes) = "lsodes"
algname(::lsodar) = "lsodar"
algname(::vode) = "vode"
algname(::daspk) = "daspk"
algname(::euler) = "euler"
algname(::rk4) = "rk4"
algname(::ode23) = "ode23"
algname(::ode45) = "ode45"
algname(::radau) = "radau"
algname(::bdf) = "bdf"
algname(::bdf_d) = "bdf_d"
algname(::adams) = "adams"
algname(::impAdams) = "impAdams"
algname(::impAdams_d) = "impAdams_d"
algname(::iteration) = "iteration"

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
    __saveat = if _saveat isa AbstractVector
        _saveat
    elseif _saveat isa Number
        collect(tspan[1]:_saveat:tspan[2])
    elseif _saveat isa Nothing
        eltype(tspan)[tspan[1], tspan[2]]
    else
        collect(_saveat)
    end

    out = rcopy(solver[].ode(times = __saveat, y = u0, func = f,
        method = algname(alg),
        parms = nothing, maxsteps = maxiters,
        rtol = reltol, atol = abstol))

    ts = @view out[:, 1]

    if u0 isa AbstractArray
        timeseries = Vector{typeof(u0)}(undef, length(ts))
        @inbounds for i in eachindex(ts)
            timeseries[i] = @view out[i, 2:end]
        end
    else
        # Scalar case: extract the solution values from column 2
        timeseries = @inbounds @view out[:, 2]
    end

    DiffEqBase.build_solution(prob, alg, ts, timeseries,
        dense = false,
        timeseries_errors = timeseries_errors)
end

@setup_workload begin
    # Precompile operations that don't require R to be initialized
    @compile_workload begin
        # Instantiate all algorithm types - this precompiles their constructors
        # and the subtype relationships
        algs = (lsoda(), lsode(), lsodes(), lsodar(), vode(), daspk(),
                euler(), rk4(), ode23(), ode45(), radau(), bdf(), bdf_d(),
                adams(), impAdams(), impAdams_d(), iteration())

        # Precompile type checks that DiffEqBase.solve will use for dispatch
        for alg in algs
            alg isa deSolveAlgorithm
            alg isa DiffEqBase.AbstractODEAlgorithm
        end

        # Precompile ODEProblem construction - common user entry point
        # This avoids ~0.3s of compilation on first ODEProblem creation
        _f_oop(u, p, t) = u
        _f_iip(du, u, p, t) = (du .= u; nothing)
        _u0 = [1.0, 0.0, 0.0]
        _tspan = (0.0, 1.0)

        # Out-of-place problem (most common)
        _prob_oop = ODEProblem(_f_oop, _u0, _tspan)
        DiffEqBase.isinplace(_prob_oop)

        # In-place problem
        _prob_iip = ODEProblem(_f_iip, _u0, _tspan)
        DiffEqBase.isinplace(_prob_iip)

        # Precompile saveat processing logic
        _saveat_arr = [0.0, 0.5, 1.0]
        _saveat_num = 0.1
        isempty(_saveat_arr)
        isempty(Float64[])
        _saveat_arr isa Array
        _saveat_num isa Number
        Array(_tspan[1]:_saveat_num:_tspan[2])

        # Precompile algorithm name extraction (used in solve)
        for alg in algs
            string(typeof(alg))[15:(end - 2)]
        end

        # Precompile array operations used in solution construction
        _out = zeros(3, 4)
        @view _out[:, 1]
        _ts = @view _out[:, 1]
        _timeseries = Vector{typeof(_u0)}(undef, length(_ts))
        for i in 1:length(_ts)
            _timeseries[i] = @view _out[i, 2:end]
        end

        # Precompile similar() for in-place function wrapper
        similar(_u0)
    end
end

end # module

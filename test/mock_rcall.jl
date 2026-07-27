using deSolveDiffEq

module MockRCall

    module MockDeSolve

        function ode(; times, y, func, parms, kwargs...)
            values = Matrix{eltype(y)}(undef, length(times), length(y) + 1)
            state = copy(y)

            for (index, time) in pairs(times)
                values[index, 1] = time
                values[index, 2:end] .= state
                state .+= func(time, state, parms)
            end

            return values
        end

    end

    rimport(::String) = MockDeSolve
    rcall(::Symbol, value) = value
    rcopy(value) = value

end

deSolveDiffEq.r_adapter[] = MockRCall
deSolveDiffEq.solver[] = nothing

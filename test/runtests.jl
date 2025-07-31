using deSolveDiffEq
using Test

function lorenz(u, p, t)
    du1 = 10.0(u[2]-u[1])
    du2 = u[1]*(28.0-u[3]) - u[2]
    du3 = u[1]*u[2] - (8/3)*u[3]
    [du1, du2, du3]
end
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob, deSolveDiffEq.lsoda())
sol = solve(prob, deSolveDiffEq.lsode())
sol = solve(prob, deSolveDiffEq.lsodes())
sol = solve(prob, deSolveDiffEq.lsodar())
sol = solve(prob, deSolveDiffEq.vode())
sol = solve(prob, deSolveDiffEq.daspk())
sol = solve(prob, deSolveDiffEq.euler())
sol = solve(prob, deSolveDiffEq.rk4())
sol = solve(prob, deSolveDiffEq.ode23())
sol = solve(prob, deSolveDiffEq.ode45())
sol = solve(prob, deSolveDiffEq.radau())
sol = solve(prob, deSolveDiffEq.bdf())
sol = solve(prob, deSolveDiffEq.bdf_d())
sol = solve(prob, deSolveDiffEq.adams())
sol = solve(prob, deSolveDiffEq.impAdams())
sol = solve(prob, deSolveDiffEq.impAdams_d())
#sol = solve(prob,deSolveDiffEq.iteration())

sol = solve(prob, deSolveDiffEq.lsoda(), saveat = 0.1)
#using Plots; plot(sol,vars=(1,2,3))

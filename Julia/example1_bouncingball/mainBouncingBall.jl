push!(LOAD_PATH, "./example1_bouncingball/dynamics")
using BouncingBall
push!(LOAD_PATH, "./utilities")
using Utilities
using LinearAlgebra
using OrdinaryDiffEq
using Plots
using LaTeXStrings

# Numerical Settings
ODESettings = Dict("solver" => Tsit5(), "reltol" =>1e-7, "abstol" => 1e-8)

analytic_gradients = true

finite_difference_step_size = 1e-7
function_tolerance = 1e-6

# Parameters
params = Dict(
    :m => 1.0,   # total mass [units don't matter, could be kg]
    :g => 1.0,   # gravity [units don't matter, could be m/s^2]
    :cor => 1.0  # coefficient of restitution [unitless]
)

# uncomment to compute the (state-based) Poincaré map
#= 
# Energy level
H_bar = 1.0

# State variables
dy = 0.0
y = H_bar

# State vector
x0 = [y; dy]

# Auxiliary variables
xi = 0.0

# Poincaré map
xt, solutions = P(x0, xi, params, ODESettings);
t2 = solutions[2].t[end] - solutions[1].t[end]
t1 = solutions[1].t[end]
xTraj2 = convert(Array, solutions[2])
x2 = xTraj2[:,1]
x1 = x0
u_init = [t2; x2; t1; x1; xi; H_bar]

# plot trajectory
plot_object = plot(title="Bouncing Ball Trajectory", xlabel=L"$t$", minorgrid=true)
for (i, sol) in enumerate(solutions)
    xTraj = convert(Array, sol)
    if i < length(solutions)
        plot!(plot_object, sol.t, xTraj[1, :], lc=:blue, label=nothing)
        plot!(plot_object, sol.t, xTraj[2, :], lc=:black, label=nothing)
    else
        plot!(plot_object, sol.t, xTraj[1, :], lc=:blue, label=L"$y$")
        plot!(plot_object, sol.t, xTraj[2, :], lc=:black, label=L"$\dot{y}$")
    end
end
display(plot_object) =#

u_init = zeros(8,1)

# Set up continuation
res(u) = res_fun(u, params, ODESettings)

opts = Dict(
    :Grad1 => analytic_gradients,
    :aimOnTarget => false,
    :idxConPar => 8,
    :FiniteDifferenceStepSize => finite_difference_step_size,
    :FunctionTolerance => function_tolerance,
    :MaxIterations => 300,
    :StepTolerance => 1e-11
)

u_init, _, _, output, jac_init = newtons_method(res, u_init, opts)

HTarget = [] # if empty, opts[:MaxIterations] points are generated
h = 1e-2
d = 1

if norm(u_init) == 0
    tangVec,_ = bifurcation_equation(u_init, res, finite_difference_step_size)
    idx = findall(abs.(tangVec[1, :]) .> function_tolerance)
    u_init .= u_init .+ h * tangVec[:, idx] .* sign.(tangVec[1, idx])
    u_init, _, _, output, jac_init = newtons_method(res, u_init, opts)
end

# Generate family of orbits
# Pre-allocate continuation data
con = Dict(
:u => zeros(length(u_init), opts[:MaxIterations]),
:tang => zeros(length(u_init), opts[:MaxIterations]),
:detSim => zeros(opts[:MaxIterations],1)',
:detAug => zeros(opts[:MaxIterations],1)'
)
numerical_continuation(con, u_init, jac_init, h, d, HTarget, res, opts)

# Plots
plot_object_1 = plot(title="Indicators for bifurcations", xlabel=L"$\bar{H}$", minorgrid=true)
plot!(plot_object_1, con[:u][end, :], con[:detSim]', label="simple Jacobian")
plot!(plot_object_1, con[:u][end, :], con[:detAug]', label="augmented Jacobian")

plot_object_2 = plot(title="Family of Normal Conservative Orbits", minorgrid=true)
plot!(plot_object_2, con[:u][end, :], con[:u][1, :] .+ con[:u][4, :])
xlabel!(L"$\bar{H}$")
ylabel!(L"$T$")

# display(plot_object_1) 
display(plot_object_2) 
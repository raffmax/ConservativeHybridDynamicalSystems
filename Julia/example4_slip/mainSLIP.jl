push!(LOAD_PATH, "./example4_slip/dynamics")
using SLIP
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
    :m    => 1.0,   # total mass [units don't matter, could be kg]
    :l_0  => 1.0,   # distance to center of mass [units don't matter, could be m]
    :g    => 1.0,   # gravity [units don't matter, could be m/s^2]
    :k    => 40,    # spring stiffness of leg [units of mg/l_0]
    :w2   => 5      # squared swing frequency [g/l_0] 
)

# state variables at simple bifurcation
l      = params[:l_0]-2*params[:m]*params[:g]/params[:k]
dl     = 0
alpha  = 0
dalpha = 0
# state vector: x_Stance = [alpha l dalpha dl]
x_sb  = [alpha;l;dalpha;dl]
# energy level
H_bar,_ = H_S(x_sb,params)

x0  = [alpha;l;dalpha;dl]

# auxiliary variables
xi  = 0; # design parameter

T = 2*pi/sqrt(params[:k]/params[:m])
t1 = T/2
t2 = 0
t3 = T/2
x1 = x0
x2 = [params[:l_0];0;0;0;0]
x3 = delta_td(x2, params) # x2->x3 in 0 time

u_init = [t3; x3; t2; x2; t1; x1; xi; H_bar]

# Set up continuation
res(u) = res_fun(u, params, ODESettings)

opts = Dict(
    :Grad1 => analytic_gradients,
    :aimOnTarget => false,
    :idxConPar => 18,
    :FiniteDifferenceStepSize => finite_difference_step_size,
    :FunctionTolerance => function_tolerance,
    :MaxIterations => 400,
    :StepTolerance => 1e-11
)

u_init, _, _, output, jac_init = newtons_method(res, u_init, opts)

HTarget = 2
h = 5e-2
d = 1

if u_init[3]==1 && u_init[6]==0
    # apply Lyapunov-Schmidt reduction to branch off into a normal orbit
    tangVec,_ = bifurcation_equation(u_init,res,finite_difference_step_size)
    idx = findall(abs.(tangVec[6,:]).>function_tolerance) # get direction into normal orbit
    u_init .= u_init .+ h * tangVec[:, idx] .* sign.(tangVec[6, idx])
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
plot!(plot_object_2, con[:u][end, :], con[:u][1, :] .+ con[:u][6, :] .+ con[:u][12, :])
xlabel!(L"$\bar{H}$")
ylabel!(L"$T$")

# display(plot_object_1) 
display(plot_object_2)
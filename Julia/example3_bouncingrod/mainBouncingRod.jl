push!(LOAD_PATH, "./example3_bouncingrod/dynamics")
using BouncingRod
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
    :r    => 0.01   # radius of gyration [units of mg/l_0]
)

# state vector: x = [y phi dy dphi]
# u_init = [t3; x3; t2; x2; t1; x1; xi; H_bar]

u_init = zeros(17,1)

# Set up continuation
res(u) = res_fun(u, params, ODESettings)

opts = Dict(
    :Grad1 => analytic_gradients,
    :aimOnTarget => false,
    :idxConPar => 17,
    :FiniteDifferenceStepSize => finite_difference_step_size,
    :FunctionTolerance => function_tolerance,
    :MaxIterations => 10000,
    :StepTolerance => 1e-11
)

u_init, _, _, output, jac_init = newtons_method(res, u_init, opts)

HTarget = 0.6
h = 2e-2
d = 1

if norm(u_init) == 0 # check if u_init is the equilibrium
   _,_,nullVec = svd(jac_init, full = true)
   nullVec = nullVec[:,15:17]
   idx = findall(abs.(nullVec[16, :]).<function_tolerance) # remove xi direction
   nullVec = nullVec[:,idx]
   # exploit symmetry: 2*t1 = t2
   coeff_sym = [hcat(2*nullVec[1,1]-nullVec[6,1], 2*nullVec[1,2]-nullVec[6,2])
                hcat(nullVec[1,1], nullVec[1,2])]\[0;1]
   tangVec  = nullVec*coeff_sym
   tangVec .= tangVec./norm(tangVec)
   u_init  .= u_init .+ 0.1 * tangVec
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
plot!(plot_object_1, con[:u][end, :], con[:detSim]'.*1000, label="simple Jacobian")
plot!(plot_object_1, con[:u][end, :], con[:detAug]', label="augmented Jacobian")

plot_object_2 = plot(title="Family of Normal Conservative Orbits", minorgrid=true)
plot!(plot_object_2, con[:u][end, :], con[:u][1, :] .+ con[:u][6, :] .+ con[:u][11, :])
xlabel!(L"$\bar{H}$")
ylabel!(L"$T$")

# display(plot_object_1) 
display(plot_object_2)
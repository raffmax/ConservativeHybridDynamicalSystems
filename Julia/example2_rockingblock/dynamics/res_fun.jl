using RockingBlock
using OrdinaryDiffEq
using LinearAlgebra
using ForwardDiff

function res_fun(u, params, ODESettings)
    n = 2

    # Extract variables from u
    t2, x2_bar, t1, x1_bar, xi, H_bar = u[1], u[2:3], u[4], u[5:6], u[7], u[8]

    # Simulate phase flows
    f1(dx, x, xi, t) = f(dx, x, xi, t, params) # (phase 1)
    x0 = x1_bar
    if t1 == 0
        x1   = x0
        Phi1 = I(n)
        psi1 = zeros(n,1)
    else
        prob_1 = ODEProblem(f1, x0, (0.0, t1), xi)
        sol_1  = solve(prob_1, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)
        function f1_fd(z)
            _prob_1 = remake(prob_1, u0 = z[1:2], p = z[3])
            solve(_prob_1, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)[1:n, 2]
        end
        dz = ForwardDiff.jacobian(f1_fd, [x0; xi])
        x1   = sol_1(t1)
        Phi1 = dz[:,1:n]
        psi1 = dz[:,n+1]
    end
    f1_val = zeros(n,1)
    f1(f1_val, x1, xi, [])
 

    f2(dx, x, xi, t) = f(dx, x, xi, t, params) # (phase 2)
    x0 = x2_bar
    if t2 == 0
        x2   = x0
        Phi2 = I(n)
        psi2 = zeros(n,1)
    else
        prob_2 = ODEProblem(f2, x0, (0.0, t2), xi)
        sol_2 = solve(prob_2, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)
        function f2_fd(z)
            _prob_2 = remake(prob_2, u0 = z[1:2], p = z[3])
            solve(_prob_2, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)[1:n, 2]
        end
        dz = ForwardDiff.jacobian(f2_fd, [x0; xi])
        x2   = sol_2(t2)
        Phi2 = dz[:,1:n]
        psi2 = dz[:,n+1]
    end
    f2_val = zeros(n,1)
    f2(f2_val, x2, xi, [])

    # Compute reset maps
    delta_12_(x) = delta_12(x, params)
    x1_P = delta_12_(x1)
    D12  = ForwardDiff.jacobian(delta_12_, x1)

    # Compute events and anchor
    e_12_(x) = e_12(x, [], [], params)
    e12_val  = e_12_(x1)
    de12 = ForwardDiff.gradient(e_12_, x1)'
    anchor_(x) = anchor(x, [], [], params)
    a  = anchor_(x2)
    da = ForwardDiff.gradient(anchor_, x2)'

    # Compute first integral
    H1, gradH1 = H(x1_bar, params)

    r = [x2 - x1_bar;
         a;
         x1_P - x2_bar;
         e12_val;
         H1 - H_bar]

    # Compute Jacobian
    R2 = [hcat(f2_val, Phi2) ; hcat(da*f2_val, da*Phi2)]
    R1 = [hcat(D12*f1_val, D12*Phi1); hcat(de12*f1_val, de12*Phi1)]
    Xi2 = vcat(psi2, da*psi2)
    Xi1 = vcat(D12*psi1, de12*psi1)
    C = [hcat(zeros(n, 1), -I(n)); zeros(1, n+1)]
    #Om = zeros(n+1, n+1)
    Or = zeros(1, n+1)
    Oc = zeros(n+1, 1)
    h1 = hcat(0, gradH1')

    dr = [hcat(R2, C, Xi2, Oc); 
          hcat(C, R1, Xi1, Oc);
          hcat(Or, h1, 0, -1)]

    return r, dr
end

using SLIP
using OrdinaryDiffEq
using LinearAlgebra
using ForwardDiff

function res_fun(u, params, ODESettings)
    nS = 4
    nF = 5

    # Extract variables from u
    t3, x3_bar, t2, x2_bar, t1, x1_bar, xi, H_bar = u[1], u[2:5], u[6], u[7:11], u[12], u[13:16], u[17], u[18]

    # Simulate phase flows
    n = nS
    f1(dx, x, xi, t) = f_S(dx, x, xi, t, params) # (phase 1)
    x0 = x1_bar
    if t1 == 0
        x1   = x0
        Phi1 = I(n)
        psi1 = zeros(n,1)
    else
        prob_1 = ODEProblem(f1, x0, (0.0, t1), xi)
        sol_1  = solve(prob_1, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)
        function f1_fd(z)
            _prob_1 = remake(prob_1, u0 = z[1:n], p = z[n+1])
            solve(_prob_1, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)[1:n, 2]
        end
        dz = ForwardDiff.jacobian(f1_fd, [x0; xi])
        x1   = sol_1(t1)
        Phi1 = dz[:,1:n]
        psi1 = dz[:,n+1]
    end
    f1_val = zeros(n,1)
    f1(f1_val, x1, xi, [])
 
    n = nF
    f2(dx, x, xi, t) = f_F(dx, x, xi, t, params) # (phase 2)
    x0 = x2_bar
    if t2 == 0
        x2   = x0
        Phi2 = I(n)
        psi2 = zeros(n,1)
    else
        prob_2 = ODEProblem(f2, x0, (0.0, t2), xi)
        sol_2 = solve(prob_2, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)
        function f2_fd(z)
            _prob_2 = remake(prob_2, u0 = z[1:n], p = z[n+1])
            solve(_prob_2, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)[1:n, 2]
        end
        dz = ForwardDiff.jacobian(f2_fd, [x0; xi])
        x2   = sol_2(t2)
        Phi2 = dz[:,1:n]
        psi2 = dz[:,n+1]
    end
    f2_val = zeros(n,1)
    f2(f2_val, x2, xi, [])

    n = nS
    f3(dx, x, xi, t) = f_S(dx, x, xi, t, params) # (phase 3)
    x0 = x3_bar
    if t3 == 0
        x3   = x0
        Phi3 = I(n)
        psi3 = zeros(n,1)
    else
        prob_3 = ODEProblem(f3, x0, (0.0, t3), xi)
        sol_3 = solve(prob_3, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)
        function f3_fd(z)
            _prob_3 = remake(prob_3, u0 = z[1:n], p = z[n+1])
            solve(_prob_3, ODESettings["solver"], reltol=ODESettings["reltol"], abstol=ODESettings["abstol"], save_everystep=false)[1:n, 2]
        end
        dz = ForwardDiff.jacobian(f3_fd, [x0; xi])
        x3   = sol_3(t3)
        Phi3 = dz[:,1:n]
        psi3 = dz[:,n+1]
    end
    f3_val = zeros(n,1)
    f3(f3_val, x3, xi, [])

    # Compute reset maps
    delta_12_(x) = delta_lo(x, params)
    x1_P = delta_12_(x1)
    D12  = ForwardDiff.jacobian(delta_12_, x1)
    delta_23_(x) = delta_td(x, params)
    x2_P = delta_23_(x2)
    D23  = ForwardDiff.jacobian(delta_23_, x2)

    # Compute events and anchor
    e_12_(x) = e_lo(x, [], [], params)
    e12_val  = e_12_(x1)
    de12 = ForwardDiff.gradient(e_12_, x1)'
    e_23_(x) = e_td(x, [], [], params)
    e23_val  = e_23_(x2)
    de23 = ForwardDiff.gradient(e_23_, x2)'
    anchor_(x) = anchor(x, [], [], params)
    a  = anchor_(x3)
    da = ForwardDiff.gradient(anchor_, x3)'

    # Compute first integral
    H1, gradH1 = H_S(x1_bar, params)

    r  = [x3-x1_bar;
          a;
          x2_P-x3_bar;
          e23_val;
          x1_P-x2_bar;
          e12_val;
          H1-H_bar]

    # Compute Jacobian
    R3 = [hcat(f3_val, Phi3) ; hcat(da*f3_val, da*Phi3)]
    R2 = [hcat(D23*f2_val, D23*Phi2) ; hcat(de23*f2_val, de23*Phi2)]
    R1 = [hcat(D12*f1_val, D12*Phi1); hcat(de12*f1_val, de12*Phi1)]
    Xi3 = vcat(psi3, da*psi3)
    Xi2 = vcat(D23*psi2, de23*psi2)
    Xi1 = vcat(D12*psi1, de12*psi1)

    CF  = [hcat(zeros(nF,1),-I(nF)); zeros(1,nF+1)]
    CS  = [hcat(zeros(nS,1),-I(nS)); zeros(1,nS+1)]
    OSF = zeros(nS+1,nF+1)   # zero matrix
    OFS = zeros(nF+1,nS+1)   # zero matrix
    OSS = zeros(nS+1,nS+1)   # zero matrix
    #OFF = zeros(nF+1,nF+1)   # zero matrix
    OcF  = zeros(nF+1,1) # zero column
    OcS  = zeros(nS+1,1) # zero column
    OrF  = zeros(1,nF+1) # zero row
    OrS  = zeros(1,nS+1) # zero row
    h1   = hcat(0, gradH1')

    dr = [hcat(R3,  OSF, CS , Xi3, OcS); 
          hcat(CS,  R2,  OSS, Xi2, OcS); 
          hcat(OFS, CF , R1,  Xi1, OcF);
          hcat(OrS, OrF, h1,  0,   -1)]

    return r, dr
end
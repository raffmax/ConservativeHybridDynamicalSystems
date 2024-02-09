using SLIP

function f_S(dx, x, xi, t, params)
    # F_S Computes the flow map of the slip 
    # during stance (i.e., phase 1 and 3)
    #
    # XDOT = F_S(X, XI, PARAMS) returns the velocity and acceleration 
    # at state X and with model parameters PARAMS.
  
    # x_Stance = [alpha l dalpha dl]

    # Parameters
    l_0  = params[:l_0]
    m    = params[:m]
    g    = params[:g]
    k    = params[:k]

    _, gradH = H_S(x, params)
    dx .= [x[3:4]; (-2*x[3]*x[4]+g*sin(x[1]))/x[2]; x[2]*x[3]^2-g*cos(x[1])-k/m*(x[2]-l_0)] + xi * gradH
end
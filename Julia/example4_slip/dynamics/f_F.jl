using SLIP

function f_F(dx, x, xi, t, params)
    # F_F Computes the flow map of the slip 
    # during flight (i.e., phase 2)
    #
    # XDOT = F_F(X, XI, PARAMS) returns the velocity and acceleration 
    # at state X and with model parameters PARAMS.
  
    # x_Flight = [y alpha dx dy dalpha]

    # Parameters
    g    = params[:g]
    w2   = params[:w2]

    _, gradH = H_F(x, params)
    dx .= [x[4:5]; 0; -g; -w2*x[2]] + xi * gradH
end
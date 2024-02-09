using RockingBlock

function f(dx, x, xi, t, params)
    # F Computes the flow map of the rocking block
    # during flight (i.e., phase 1)
    #
    # F(DX, X, XI, T, PARAMS) returns the velocity and acceleration 
    # at state X and with model parameters PARAMS.

    p     = 3/4*params[:g]/params[:l_0];
    beta  = params[:beta];

    _, gradH = H(x, params)
    dx .= [x[2]; -p*sin(beta-x[1])] + xi * gradH
end
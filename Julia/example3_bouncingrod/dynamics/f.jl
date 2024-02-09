using BouncingRod

function f(dx, x, xi, t, params)
    # F Computes the flow map of the bouncing rod
    # during flight (i.e., phase 1, 2, 3)
    #
    # F(DX, X, XI, T, PARAMS) returns the velocity and acceleration 
    # at state X and with model parameters PARAMS.

    _, gradH = H(x, params)
    dx .= [x[3:4];-params[:g]; 0] + xi * gradH
end
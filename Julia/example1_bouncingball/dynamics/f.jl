using BouncingBall

function f(dx, x, xi, t, params)
    # F Computes the flow map of the bouncing ball
    # during flight (i.e., phase 1)
    #
    # F(DX, X, XI, T, PARAMS) returns the velocity and acceleration 
    # at state X and with model parameters PARAMS.

    ~, gradH = H(x, params)
    dx .= [x[2]; -params[:g]] + xi * gradH
end
function e_23(x, t, integrator, params)
    # E_23 Computes the right touch-down event for the bouncing rod
    #
    # V = E_23(X,PARAMS) returns condition for event 
    # V is zero at touchdown.

    # Event is detected if touch-down
    v = x[1]+params[:l_0]*sin(x[2])/2

    return v
end
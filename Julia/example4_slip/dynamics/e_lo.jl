function e_lo(x, t, integrator, params)
    # E_LO Computes the lift-off event for the slip model
    #
    # V = E_LO(X,PARAMS) returns condition for event 
    # V is zero at touchdown.

    # Event is detected if touch-down
    v = params[:l_0]-x[2]

    return v
end
function e_12(x, t, integrator, params)
    # E_12 Computes the jump set for the rocking block
    #
    # V = E_12(X,PARAMS) returns condition for event 
    # V is zero at touchdown.

    # Event is detected if touch-down
    v = x[1]

    return v
end
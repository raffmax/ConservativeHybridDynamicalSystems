function e_td(x, t, integrator, params)
    # E_TD Computes the jump set for the bouncing ball
    #
    # V = E_TD(X,PARAMS) returns condition for event 
    # V is zero at touchdown.

    # Event is detected if touch-down
    v = x[1]

    return v
end
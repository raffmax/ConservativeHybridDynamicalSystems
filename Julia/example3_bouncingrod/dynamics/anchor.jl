function anchor(x, t, integrator, params)
    # A Computes the anchor for the bouncing rod
    #
    # V = ANCHOR(X,PARAMS) returns V, ISTERM, and DIR 
    # according to the ode45 event function documentation. V is zero at apex.

    # Event is detected if apex
    v = x[3]

    return v
end
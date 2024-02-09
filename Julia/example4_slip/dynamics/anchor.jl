function anchor(x, t, integrator, params)
    # A Computes the anchor for the slip model
    #
    # V = ANCHOR(X,PARAMS) returns V, ISTERM, and DIR 
    # according to the ode45 event function documentation. V is zero at nadir.

    # Event is detected if nadir
    v = -x[4]

    return v
end
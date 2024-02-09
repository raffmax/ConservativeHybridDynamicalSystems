function e_td(x, t, integrator, params)
    # E_TD Computes the touch-down event for the slip model
    #
    # V = E_TD(X,PARAMS) returns condition for event 
    # V is zero at touchdown.

    # Event is detected if touch-down
    v = x[1]-params[:l_0]*cos(x[2])

    return v
end
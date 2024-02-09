using BouncingBall
using OrdinaryDiffEq
affect!(integrator) = terminate!(integrator) # terminate all ODEs after detecting an event

function P(xk, xi, params, ODESettings)
    # P Compute the stride-to-stride map of the BOUNCING BALL.
    #   XKPLUS1 = P(XK, PARAMS) returns the state XKPLUS1 after 
    #   a single step that started from XK based on the flow map, 
    #   jump map, and jump set, for a model with parameters PARAMS.
    #   The phase sequence is 1-2:
    #   1: flight post apex
    #   2: flight pre apex

    # start time
    t0 = 0.0

    # max simulation time to wait for all events to happen for a complete
    # stride
    t_max = 50.0

    ## Compute a stride
    # A stride consists of successful transition from a fixed contact sequence:
    # phase 1 -> phase 1
    #  xk -------o---------o-------> xkPLUS1

    ## Compute flow in phase 1 (flight post apex)
    ode_fun_1(dx, x, xi, t) = f(dx, x, xi, t, params)  # flight (phase 1)
    ode_event_1(x, t, integrator) =  e_td(x, t, integrator, params)   # Touch Down
    cb_1 = ContinuousCallback(ode_event_1, affect!) 
    prob_1 = ODEProblem(ode_fun_1,xk,(t0,t0+t_max),xi)
    sol_1 = solve(prob_1, ODESettings["solver"], saveat=0.01, callback=cb_1, reltol=ODESettings["reltol"], abstol=ODESettings["abstol"])

    xe_prek = sol_1[end]
    tek     = sol_1.t

    # apply jump map?
    if tek[end] < t0+t_max
        t0 = tek[end]
        
        # jump map
        xk = delta_td(xe_prek, params)
        
    else
        warn("No collision for $t0 <= t <= $t_max.")
        # return empty list to signal that we have not completed a stride
        return [], []
    end
    
    ## Compute flow in phase 2 (flight pre apex)
    ode_fun_2(dx, x, xi, t) = f(dx, x, xi, t, params)  # flight (phase 2)
    ode_event_2(x, t, integrator) =  anchor(x, t, integrator, params)   # Apex Transit
    cb_2 = ContinuousCallback(ode_event_2, affect!)
    prob_2 = ODEProblem(ode_fun_2,xk,(t0,t0+t_max),xi)
    sol_2 = solve(prob_2, ODESettings["solver"], saveat=0.01, callback=cb_2, reltol=ODESettings["reltol"], abstol=ODESettings["abstol"])

    xe_prek = sol_2[end]
    xkPLUS1 = xe_prek

    # create trajectory struct for output
    solutions = [sol_1,sol_2]
    return xkPLUS1, solutions
end

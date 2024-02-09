function numerical_continuation(con, u_init, jac_init, h, d, HTarget, resFcn, opts)
    tang_init = get_tangent(jac_init)

    # Store initial point in continuation data
    con[:u][:, 1]    .= u_init
    con[:tang][:, 1] .= tang_init
    con[:detSim][1]   = det(jac_init[:, 1:end-1])
    con[:detAug][1]   = det([jac_init; tang_init'])

    # initialize continuation
    u = copy(u_init)
    tang = copy(tang_init)

    if !isempty(HTarget)
        dirCon = sign(HTarget - u_init[end])
        # define correct direction of continuation
        if sign(HTarget - u_init[end]) * sign(tang_init[end]) == -1
            d = -1
        end
    else
        HTarget = 0
        dirCon  = 0
    end

    endLoop = false
    iter = 1

    while !endLoop
        # predictor
        u_pred = u + d * h * tang
        if dirCon * u_pred[end] > dirCon * HTarget
            endLoop = true
            opts[:aimOnTarget] = true
            u_pred .= u_pred .- (u_pred[end] - HTarget) .* tang ./ tang[end]
        end

        # corrector
        u_corr, fval, exitflag, output, jac = newtons_method(resFcn, u_pred, opts)
        tang_new = get_tangent(jac)

        # check for bifurcation
        if dot(tang_new, tang) < 0
            println("Simple Bifurcation between H=", u[end], " and H=", u_pred[end])
            d = -d
        elseif tang_new[end]*tang[end] < 0
            println("Turning Point between H=", u[end], " and H=", u_pred[end])
        end

        # update
        u = copy(u_corr)
        tang = copy(tang_new)

        if dirCon * u_corr[end] > dirCon * HTarget
            # the corrector found a solution with H > HTarget
            # however, the predictor was H < HTarget
        else
            iter += 1
            # Store continuation data
            con[:u][:, iter]    .= u
            con[:tang][:, iter] .= tang
            con[:detSim][iter]   = det(jac[:, 1:end-1])
            con[:detAug][iter]   = det([jac; tang'])
        end

        if endLoop
            # iter <= opts[:MaxIterations]
            con[:u]      = con[:u][:, 1:iter]
            con[:tang]   = con[:tang][:, 1:iter]
            con[:detSim] = con[:detSim][:, 1:iter]
            con[:detAug] = con[:detAug][:, 1:iter]
        end

        # terminate when max iterations
        if iter == opts[:MaxIterations]
            endLoop = true
        end
    end
end

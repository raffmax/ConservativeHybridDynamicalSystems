function get_jacobian_fd(fcn, x, step_size)
    # get f(x)
    f = fcn(x)
    
    # get dimensions
    nF = length(f)
    nX = length(x)

    jacobian = zeros(nF, nX)

    for iCol in 1:nX
        iX = zeros(nX)
        iX[iCol] = 1
        f_h = fcn(x + h * iX)
        jacobian[:, iCol] = (f_h - f) / step_size
    end

    return f, jacobian
end
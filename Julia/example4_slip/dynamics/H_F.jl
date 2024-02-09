function H_F(x, params)
    # H_F First integral of slip model during flight

    # Parameters
    m    = params[:m]
    g    = params[:g]

    # x_Flight = [y alpha dx dy dalpha]

    # Calculate value of first integral
    v = m*(g*x[1]+0.5*(x[3]^2+x[4]^2))
    
    # Calculate first derivative of first integral: dH
    gradH = m*[g;0;x[3];x[4];0]
    
    return v, gradH
end
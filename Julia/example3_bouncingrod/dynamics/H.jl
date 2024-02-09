function H(x, params)
    # Access parameters from the dictionary

    r    = params[:r]
    m    = params[:m]
    g    = params[:g]

    # Calculate value of first integral
    v = m*g*x[1]+m*(x[3]^2+r^2*x[4]^2)/2
    
    # Calculate first derivative of first integral: dH
    gradH = m*vcat(g, 0, x[3], r^2*x[4])
    
    return v, gradH
end
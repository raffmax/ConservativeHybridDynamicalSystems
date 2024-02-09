function H(x, params)
    # Access parameters from the dictionary
    m = params[:m]
    g = params[:g]

    # Calculate value of first integral
    v = m * (g * x[1] + 0.5 * x[2]^2)
    
    # Calculate first derivative of first integral: dH
    gradH = m * [g; x[2]]
    
    return v, gradH
end
function H(x, params)
    # Access parameters from the dictionary

    R    = params[:l_0]
    beta = params[:beta]
    m    = params[:m]
    g    = params[:g]

    # Calculate value of first integral
    v = 2/3*m*R^2*x[2]^2+R*(cos(x[1]-beta)-cos(beta))*m*g
    
    # Calculate first derivative of first integral: dH
    gradH = [-R*m*g*sin(x[1]-beta) ; 4/3*m*R^2*x[2]]
    
    return v, gradH
end
function H_S(x, params)
    # H_S First integral of slip model during stance

    # Parameters
    l_0  = params[:l_0]
    m    = params[:m]
    g    = params[:g]
    k    = params[:k]

    # x_Stance = [alpha l dalpha dl]

    # Calculate value of first integral
    v = m*g*cos(x[1])*x[2]+0.5*m*(x[2]^2*x[3]^2+x[4]^2)+0.5*k*(x[2]-l_0)^2
    
    # Calculate first derivative of first integral: dH
    gradH = m*vcat(-g*sin(x[1])*x[2], 
                    g*cos(x[1])+m*x[2]*x[3]^2+k*(x[2]-l_0),
                    x[2]^2*x[3],
                    x[4])
    
    return v, gradH
end
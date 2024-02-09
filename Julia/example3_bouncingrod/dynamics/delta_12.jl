function delta_12(x_pre, params)
    # delta_12 Computes the reset map of the bouncing rod at left touch-down.

    # DELTA_12(X_PRE, PARAMS) returns the post-impact state 
    # X_POST given pre-impact state X_PRE for model with parameters PARAMS.

    # x = [y phi dy dphi] 
    phi  = x_pre[2]
    dy   = x_pre[3]
    dphi = x_pre[4]
    
    c = cos(phi)
    s = sin(phi)
    l = params[:l_0]
    r = params[:r]
    
    vec    = vcat(0, 0, -l*dphi*r^2*c+2*dy*r^2, -l*dy*c+(l^2*dphi*c^2)/2)
    x_post = x_pre - 4*vec/(l^2*c^2+4*r^2)

    return x_post
end
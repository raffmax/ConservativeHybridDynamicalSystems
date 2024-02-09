function delta_td(x_pre, params)
    # delta_td Computes the jump map of the BOUNCING BALL.

    # DELTA_TD(X_PRE, PARAMS) returns the post-impact state 
    # X_POST given pre-impact state X_PRE for model with parameters PARAMS.

    D      = [1 0; 0 -params[:cor]]
    x_post = D * x_pre

    return x_post
end
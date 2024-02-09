function [x_post,D] = delta_td(x_pre, params)
% delta_td Computes the jump map of the BOUNCING BALL.
%
% [X_POST, D] = DELTA_TD(X_PRE, PARAMS) returns the post-impact state 
% X_POST given pre-impact state X_PRE for model with parameters PARAMS.

cor = params.cor;

D      = [1 0; 0 -cor];
x_post = D * x_pre;
end
function [x_post,D] = delta_12(x_pre, params)
% delta_12 Computes the jump map of the rocking block at left touch-down.
%
% [X_POST, D] = DELTA_12(X_PRE, PARAMS) returns the post-impact state 
% X_POST given pre-impact state X_PRE for model with parameters PARAMS.

% x = [phi dphi] 
x_post = [1,0;0,-1]*x_pre;

D = [1,0;0,-1];
end
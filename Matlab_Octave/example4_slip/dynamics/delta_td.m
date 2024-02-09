function [x_post,D] = delta_td(x_pre, params)
% delta_td Computes the jump map of the slip at touch-down.
%
% [X_POST, D] = DELTA_TD(X_PRE, PARAMS) returns the post-impact state 
% X_POST given pre-impact state X_PRE for model with parameters PARAMS.

% x_Flight = [y alpha dx dy dalpha] -> x_pre
% x_Stance = [alpha l dalpha dl] -> x_post
x = x_pre;

% Parameters
l_0 = params.l_0;

x_post = [x(2);
          l_0;
          -(cos(x(2))*x(3)+sin(x(2))*x(4))/l_0;
          cos(x(2))*x(4)-sin(x(2))*x(3)];

D = [0,1,0,0,0;...
     0,0,0,0,0;...
     0,(sin(x(2))*x(3)-cos(x(2))*x(4))/l_0,-cos(x(2))/l_0,-sin(x(2))/l_0,0;...
     0,-sin(x(2))*x(4)-cos(x(2))*x(3),-sin(x(2)),cos(x(2)),0];
end
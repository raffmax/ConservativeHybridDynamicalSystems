function [x_post,D] = delta_lo(x_pre, params)
% delta_lo Computes the jump map of the slip at lift-off.
%
% [X_POST, D] = DELTA_LO(X_PRE, PARAMS) returns the post-impact state 
% X_POST given pre-impact state X_PRE for model with parameters PARAMS.

% x_Stance = [alpha l dalpha dl] -> x_pre
% x_Flight = [y alpha dx dy dalpha] -> x_post
x = x_pre;

x_post = [x(2)*cos(x(1));
          x(1);
          -x(4)*sin(x(1))-x(2)*x(3)*cos(x(1));
          x(4)*cos(x(1))-x(2)*x(3)*sin(x(1));
          x(3)];

D = [-x(2)*sin(x(1)),cos(x(1)),0,0;...
     1,0,0,0;...
     -x(4)*cos(x(1))+x(2)*x(3)*sin(x(1)),-x(3)*cos(x(1)),-x(2)*cos(x(1)),-sin(x(1));...
     -x(4)*sin(x(1))-x(2)*x(3)*cos(x(1)),-x(3)*sin(x(1)),-x(2)*sin(x(1)),cos(x(1));...
     0,0,1,0];
end
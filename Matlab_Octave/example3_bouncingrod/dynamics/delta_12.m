function [x_post,D] = delta_12(x_pre, params)
% delta_12 Computes the reset map of the bouncing rod at left touch-down.
%
% [X_POST, D] = DELTA_12(X_PRE, PARAMS) returns the post-impact state 
% X_POST given pre-impact state X_PRE for model with parameters PARAMS.

% x = [y phi dy dphi] 
phi  = x_pre(2);
dy   = x_pre(3);
dphi = x_pre(4);

c = cos(phi);
s = sin(phi);
l = params.l_0;
r = params.r;

vec    = [0;
          0;
          -l*dphi*r^2*c + 2*dy*r^2;
          -l*dy*c + (l^2*dphi*c^2)/2];
x_post = x_pre - 4*vec/(l^2*c^2+4*r^2);

D = [[1,                                                                                              0,                                       0,                                               0];...
[0,                                                                                                   1,                                       0,                                               0];...
[0, -(4*l*r^2*s*(4*dphi*r^2 + 4*dy*l*c - dphi*l^2*c^2))/(l^2*c^2 + 4*r^2)^2,    1 - (8*r^2)/(l^2*c^2 + 4*r^2),     (4*l*r^2*c)/(l^2*c^2 + 4*r^2)];...
[0,    (4*l*s*(dy*l^2*c^2 - 4*dy*r^2 + 4*dphi*l*r^2*c))/(l^2*c^2 + 4*r^2)^2, (4*l*c)/(l^2*c^2 + 4*r^2), 1 - (2*l^2*c^2)/(l^2*c^2 + 4*r^2)]];
end
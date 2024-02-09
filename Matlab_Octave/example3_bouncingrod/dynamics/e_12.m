function [v, isterm, dir, de] = e_12(x, params)
% E_12 Computes the left touch-down event for the bouncing rod
%
% [V, ISTERM, DIR, DE] = E_12(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at 
% left touch down.

% Event is detected if touch-down
v = x(1)-params.l_0*sin(x(2))/2;

isterm = 1;
dir    = -1;

de     = [1 -params.l_0*cos(x(2))/2 0 0];
end
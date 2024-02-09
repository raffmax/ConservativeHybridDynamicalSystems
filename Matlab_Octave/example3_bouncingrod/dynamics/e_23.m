function [v, isterm, dir, de] = e_23(x, params)
% E_23 Computes the right touch-down event for the bouncing rod
%
% [V, ISTERM, DIR, DE] = E_23(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at 
% right touch down.

% Event is detected if touch-down
v = x(1)+params.l_0*sin(x(2))/2;

isterm = 1;
dir    = -1;

de     = [1 params.l_0*cos(x(2))/2 0 0];
end
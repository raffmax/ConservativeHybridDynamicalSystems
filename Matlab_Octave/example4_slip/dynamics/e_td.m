function [v, isterm, dir, de] = e_td(x, params)
% E_TD Computes the jump set for the slip
%
% [V, ISTERM, DIR, DE] = E_TD(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at 
% touchdown.

% x_Flight = [y alpha dx dy dalpha]
% Event is detected if touch-down
v = x(1)-params.l_0*cos(x(2));

isterm = 1;
dir    = -1;

de     = [1 params.l_0*sin(x(2)) 0 0 0];
end
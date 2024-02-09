function [v, isterm, dir, de] = e_td(x, params)
% E_TD Computes the jump set for the bouncing ball
%
% [V, ISTERM, DIR, DE] = E_TD(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at touchdown.

% Event is detected if touch-down
v = x(1);

isterm = 1;
dir    = -1;

de     = [1 0];
end
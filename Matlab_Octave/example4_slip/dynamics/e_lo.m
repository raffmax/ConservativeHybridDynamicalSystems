function [v, isterm, dir, de] = e_lo(x, params)
% E_LO Computes the jump set for the slip
%
% [V, ISTERM, DIR, DE] = E_LO(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at 
% lift-off.

% x_Stance = [alpha l dalpha dl]
% Event is detected if touch-down
v = params.l_0-x(2);

isterm = 1;
dir    = -1;

de     = [0 -1 0 0];
end
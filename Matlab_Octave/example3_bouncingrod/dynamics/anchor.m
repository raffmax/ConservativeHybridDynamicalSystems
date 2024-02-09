function [v, isterm, dir, da] = anchor(x, params)
% A Computes the anchor for the bouncing rod
%
% [V, ISTERM, DIR, DA] = ANCHOR(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at apex.

% Event is detected if apex
v = x(3);

isterm = 1;
dir    = -1;

da     = [0,0,1,0];
end
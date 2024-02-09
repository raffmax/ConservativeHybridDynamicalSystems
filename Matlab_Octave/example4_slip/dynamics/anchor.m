function [v, isterm, dir, da] = anchor(x, params)
% A Computes the anchor for the slip model
%
% [V, ISTERM, DIR, DA] = ANCHOR(X,PARAMS) returns V, ISTERM, and DIR 
% according to the ode45 event function documentation. V is zero at nadir.

% Event is detected if nadir
v = -x(4);

isterm = 1;
dir    = -1;

da     = [0,0,0,-1];
end
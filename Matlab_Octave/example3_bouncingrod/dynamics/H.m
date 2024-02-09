function [v,dH,ddH] = H(x,params)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
v   = params.m*params.g*x(1)+params.m*(x(3)^2+params.r^2*x(4)^2)/2;
dH  = params.m*[params.g,0,x(3),params.r^2*x(4)];
ddH = params.m*diag([0,0,1,params.r^2]);
end

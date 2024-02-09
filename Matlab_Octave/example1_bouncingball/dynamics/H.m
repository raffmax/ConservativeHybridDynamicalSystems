function [v,dH,ddH] = H(x,params)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
v   = params.m*(params.g*x(1)+0.5*x(2).^2);
dH  = params.m*[params.g,x(2)];
ddH = [0 0; 0 params.m];
end

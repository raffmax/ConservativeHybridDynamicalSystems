function [v,dH,ddH] = H(x,params)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

R     = params.l_0;
beta  = params.beta;
g     = params.g;
m     = params.m;
    
v   = 2/3*m*R^2*x(2)^2+R*(cos(x(1)-beta)-cos(beta))*m*g;
dH  = [-R*m*g*sin(x(1)-beta),4/3*m*R^2*x(2)];
ddH = R*m*diag([-g*cos(x(1)-beta),4/3*R]);
end

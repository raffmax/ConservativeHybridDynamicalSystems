function [H,dH,ddH] = H_F(x,params)
%H_F First integral of slip model during flight

% Parameters
g = params.g;
m = params.m;
    
% x_Flight = [y alpha dx dy dalpha]
H   = m*(g*x(1)+0.5*(x(3).^2+x(4).^2));
dH  = m*[g,0,x(3),x(4),0];
ddH = m*diag([0,0,1,1,0]);
end

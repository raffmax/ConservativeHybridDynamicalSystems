function xdot = f_S(x,xi,params)
% F_S Computes the flow map of the slip 
% during stance (i.e., phase 1 and 3)
%
% XDOT = F_S(X, XI, PARAMS) returns the velocity and acceleration 
% at state X and with model parameters PARAMS.
n = 4;
% x_Stance = [alpha l dalpha dl]

% Parameters
g = params.g;
k = params.k;
m = params.m;
l_0 = params.l_0;


if numel(x) == n
    [~,dH] = H_S(x,params);
    gradH  = dH';
    xdot   = [x(3:4);...
              (-2*x(3)*x(4)+g*sin(x(1)))/x(2);...
              x(2)*x(3)^2-g*cos(x(1))-k/m*(x(2)-l_0)]...
             + xi*gradH;
else
    X = x; % includes sensitivities
    % reshape X
    x   = X(1:n);
    Phi = reshape(X(n+(1:n^2)),[n,n]);
    psi = X(n+n^2+(1:n));

    [~,dH,ddH] = H_S(x,params);
    gradH = dH';

    xdot   = [x(3:4);...
              (-2*x(3)*x(4)+g*sin(x(1)))/x(2);...
              x(2)*x(3)^2-g*cos(x(1))-k/m*(x(2)-l_0)]...
             + xi*gradH;
    df     = [zeros(2) eye(2);
              [g*cos(x(1)),-xdot(3),-2*x(4),-2*x(3)]/x(2);
              +g*sin(x(1)),x(3)^2-k/m,2*x(2)*x(3),0];
    Phidot = (df+xi*ddH)*Phi;
    psidot = (df+xi*ddH)*psi + gradH;

    xdot = [xdot;Phidot(:);psidot]; % includes sensitivities
end

end
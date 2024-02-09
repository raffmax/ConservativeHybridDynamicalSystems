function xdot = f(x,xi,params)
% F Computes the flow map of the rocking block
% during flight (i.e., phase 1 2 3)
%
% XDOT = F(X, XI, PARAMS) returns the velocity and acceleration 
% at state X and with model parameters PARAMS.
n = 2;

p     = 3/4*params.g/params.l_0;
beta  = params.beta;

if numel(x) == n
    [~,dH] = H(x,params);
    gradH = dH';
    xdot = [x(2);-p*sin(beta-x(1))] + xi*gradH;
else
    X = x; % includes sensitivities
    % reshape X
    x   = X(1:n);
    Phi = reshape(X(n+(1:n^2)),[n,n]);
    psi = X(n+n^2+(1:n));

    [~,dH,ddH] = H(x,params);
    gradH = dH';

    xdot = [x(2);-p*sin(beta-x(1))] + xi*gradH;
    Phidot = ([0 1;p*cos(beta-x(1)) 0]+xi*ddH)*Phi;
    psidot = ([0 1;p*cos(beta-x(1)) 0]+xi*ddH)*psi + gradH;

    xdot = [xdot;Phidot(:);psidot]; % includes sensitivities
end

end

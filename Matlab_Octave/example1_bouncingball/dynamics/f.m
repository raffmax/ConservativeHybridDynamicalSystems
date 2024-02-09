function xdot = f(x,xi,params)
% F Computes the flow map of the bouncing ball
% during flight (i.e., phase 1)
%
% XDOT = F(X, XI, PARAMS) returns the velocity and acceleration 
% at state X and with model parameters PARAMS.
n = 2;

if numel(x) == n
    [~,dH] = H(x,params);
    gradH = dH';
    xdot = [x(2); -params.g] + xi*gradH;
else
    X = x; % includes sensitivities
    % reshape X
    x   = X(1:n);
    Phi = reshape(X(n+(1:n^2)),[n,n]);
    psi = X(n+n^2+(1:n));

    [~,dH,ddH] = H(x,params);
    gradH = dH';

    xdot   = [x(2); -params.g] + xi*gradH;
    Phidot = ([0 1; 0 0]+xi*ddH)*Phi;
    psidot = ([0 1; 0 0]+xi*ddH)*psi + gradH;

    xdot = [xdot;Phidot(:);psidot]; % includes sensitivities
end

end

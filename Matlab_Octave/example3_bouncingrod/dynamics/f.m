function xdot = f(x,xi,params)
% F Computes the flow map of the rocking block
% during flight (i.e., phase 1 2 3)
%
% XDOT = F(X, XI, PARAMS) returns the velocity and acceleration 
% at state X and with model parameters PARAMS.
n = 4;

if numel(x) == n
    [~,dH] = H(x,params);
    gradH = dH';
    xdot = [x(3:4);-params.g;0] + xi*gradH;
else
    X = x; % includes sensitivities
    % reshape X
    x   = X(1:n);
    Phi = reshape(X(n+(1:n^2)),[n,n]);
    psi = X(n+n^2+(1:n));

    [~,dH,ddH] = H(x,params);
    gradH = dH';

    xdot   = [x(3:4);-params.g;0] + xi*gradH;
    Phidot = ([zeros(2),eye(2);zeros(2,4)]+xi*ddH)*Phi;
    psidot = ([zeros(2),eye(2);zeros(2,4)]+xi*ddH)*psi + gradH;

    xdot = [xdot;Phidot(:);psidot]; % includes sensitivities
end

end

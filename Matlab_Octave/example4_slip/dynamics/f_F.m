function xdot = f_F(x,xi,params)
% F_F Computes the flow map of the slip 
% during flight (i.e., phase 2)
%
% XDOT = F_F(X, XI, PARAMS) returns the velocity and acceleration 
% at state X and with model parameters PARAMS.
n = 5;
% x_Flight = [y alpha dx dy dalpha]

% Parameters
g = params.g;
w2 = params.w2;

if numel(x) == n
    [~,dH] = H_F(x,params);
    gradH  = dH';
    xdot   = [x(4:5); 0; -g; -w2*x(2)] + xi*gradH;
else
    X = x; % includes sensitivities
    % reshape X
    x   = X(1:n);
    Phi = reshape(X(n+(1:n^2)),[n,n]);
    psi = X(n+n^2+(1:n));

    [~,dH,ddH] = H_F(x,params);
    gradH = dH';

    xdot   = [x(4:5); 0; -g; -w2*x(2)] + xi*gradH;
    df     = [zeros(2,3) eye(2); zeros(2,5); 0 -w2 0 0 0];
    Phidot = (df+xi*ddH)*Phi;
    psidot = (df+xi*ddH)*psi + gradH;

    xdot = [xdot;Phidot(:);psidot]; % includes sensitivities
end

end

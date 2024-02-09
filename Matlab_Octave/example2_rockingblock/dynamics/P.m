function [xkPLUS1,traj] = P(xk, xi, params, ode)
% P Computes the Poincaré map of the rocking block.
%   XKPLUS1 = P(XK, PARAMS) returns the state XKPLUS1 after 
%   a single step that started from XK based on the flow map, 
%   jump map, and jump set, for a model with parameters PARAMS.
%   The phase sequence is 1-2:
%   1: flight post apex
%   2: flight after left and before right contact
%   3: stance pre apex
%
% See also: ODE45

%  Institute for Nonlinear Mechanics 04/12/2023, Matlab R2022a, v1

% start time
t0 = 0;

% max simulation time to wait for all events to happen for a complete
% stride
t_max = 50;

% initialize data logging of trajectory and events
tout = t0;  % log flow
xout = xk;

te      = [];    % log jumps
xe_pre  = {};
xe_post = {};
ie      = [];

% set up functions for ode45
%   we work with autonomous systems f(x), e(x), etc., but ode45 requires 
%   fun = f(t, x) and event = e(t, x), so we'll wrap our functions
%   accordingly
% since we implement 3 phases for the slip model, we have 3 individual f, g
% and e functions

%% Compute a stride
% A stride consists of successful transition from a fixed contact sequence:
% phase 1 -> phase 2 -> phase 3 -> phase 1
%  xk -------o---------o-------> xkPLUS1

%% Compute flow in phase 1 (post apex)
ode_fun_1   = @(t, x) f(x, xi, params);  % flight (phase 1)
ode_event_1 = @(t, x) e_12(x, params);   % Touch Down
opts_1 = odeset('Events', ode_event_1, 'RelTol', ode.settings.RelTol,...
              'AbsTol', ode.settings.AbsTol);
[tk, xk, tek, xe_prek, iek] = ode.solver(ode_fun_1, [t0, t_max], xk, opts_1);

% log data
% convert ode45's output into our expected output format
tout   = [tout, tk'];
xout   = [xout, xk'];
te     = [te, tek'];
xe_pre = [xe_pre, xe_prek'];
ie     = [ie, iek];

% apply jump map?
if(~isempty(tek))
    t0 = tek(end);
    
    % jump map
    xk = delta_12(xe_prek', params);
    
    % log data
    xe_post = [xe_post, xk];
else
    warning('No collision for %.2f <= t <= %.2f.', t0, t_max);
    % return empty list to signal that we have not completed a stride
    xkPLUS1 = [];
end
    
%% Compute flow in phase 2 (pre apex)
ode_fun_2   = @(t, x) f(x, xi, params); % flight (mode 3)
ode_event_2 = @(t, x) anchor(x, params); % Apex Transit
opts_2 = odeset('Events', ode_event_2, 'RelTol', ode.settings.RelTol,...
                'AbsTol', ode.settings.AbsTol);
[tk, xk, tek, xe_prek, iek] = ode45(ode_fun_2, [t0, t_max], xk, opts_2);

% log data
% convert ode45's output into our expected output format
tout   = [tout, tk'];
xout   = [xout, xk'];
te     = [te, tek'];
xe_pre = [xe_pre, xe_prek'];
ie     = [ie, iek];

% apply jump map?
if(~isempty(tek))
    t0      = tek;
    xk      = xe_prek';
    xe_post = [xe_post, xk];
    % return final post-impact state
    xkPLUS1 = xe_post(:, end);
else
    warning('No collision for %.2f <= t <= %.2f.', t0, t_max);
    % return empty list to signal that we have not completed a stride
    xkPLUS1 = [];
end

% create trajectory struct for output
traj.t       = tout;          % Row vector of the time of the trajectory
traj.x       = xout;          % Row of states (as column vectors) for each 
                              % time step
traj.te      = te;            % Row vector of the times of all events
traj.ie      = ie;            % Row vector of the ID-numbers of all events
traj.xe_pre  = xe_pre;        % Row of pre-jump states (as column vectors) 
                              % for each event
traj.xe_post = xe_post;       % Row of post-jump states (as column vectors) 
                              % for each event


end
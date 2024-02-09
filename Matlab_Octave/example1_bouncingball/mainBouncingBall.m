restoredefaultpath
addpath('./dynamics');
parentDir = fileparts(pwd);
addpath([parentDir,'/utilities']);


%% Numerical Settings
ode.settings = odeset('RelTol',1e-7,'AbsTol',1e-8);
ode.solver   = @ode45;

analyticGradients = true; % use analytic derivatives / sensitivities

FiniteDifferenceStepSize = 1e-7; % if analyticGradients == false
FunctionTolerance        = 1e-6; % tolerance of r(.)=0

%% Parameters
params       = struct();
params.m     = 1;   % total mass [units don't matter, could be kg]
params.g     = 1;   % gravity [units don't matter, could be m/s^2]
params.cor   = 1;  % coefficient of restitution [unitless]

% uncomment to compute the (state-based) Poincaré map
% 
% % energy level
% H_bar = 1;
% 
% % state variables
% dy =  0;
% y  =  H_bar;
% % state vector
% x0  = [y;dy];
% 
% % auxiliary variables
% 
% xi  = 0; % design parameter
% % Poincaré map
% [xt, traj] = P(x0, xi, params, ode);
% 
% % z_init = [t2; x2; t1; x1; xi]; % Initial conditions
% t2 = traj.te(2)-traj.te(1);
% t1 = traj.te(1);
% x2 = traj.xe_post(:,1);
% x1 = x0;
%
% u_init = [t2; x2; t1; x1; xi; H_bar];

u_init = zeros(8,1);

%% set up continuation
res = @(u) resFun(u,params,ode);

opts             = [];
opts.Grad1       = analyticGradients;
opts.aimOnTarget = false;
opts.idxConPar   = 8;
opts.FiniteDifferenceStepSize = FiniteDifferenceStepSize;
opts.FunctionTolerance        = FunctionTolerance;
opts.MaxIterations     = 300;
opts.StepTolerance     = 1e-11;
[u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);

HTarget = [];
h = 1e-2; % fixed step-size
d = +1; % direction of arclength continuation

if norm(u_init)==0 % check if u_init is the equilibrium
   % apply Lyapunov-Schmidt reduction to branch off into a normal orbit
   tangVec = BifurcationEquation(u_init,res,FiniteDifferenceStepSize);
   idx = find(abs(tangVec(1,:))>FunctionTolerance); % get direction into normal orbit
   u_init = u_init + h*tangVec(:,idx)*sign(tangVec(1,idx));
   [u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);
end

% generate family of orbits
con = numericalContinuation(u_init,jac_init,h,d,HTarget,res,opts);

%% plots
figure
grid on
hold on
plot(con.u(end,:),con.detSim)
plot(con.u(end,:),con.detAug)
xlabel('$\bar{H}$','Interpreter','latex')
legend('simple Jacobian','augmented Jacobian','Interpreter','latex')

figure
grid on
hold on
plot(con.u(end,:),con.u(1,:)+con.u(4,:))
xlabel('$\bar{H}$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')

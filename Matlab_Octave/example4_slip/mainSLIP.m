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
params.l_0   = 1;   % leg length [units don't matter, could be m]
params.g     = 1;   % gravity [units don't matter, could be m/s^2]
params.k     = 40;  % spring stiffness of leg [units of mg/l_0]
params.w2    = 5;   % squared swing frequency [g/l_0] 

% state variables at simple bifurcation
l      = params.l_0-2*params.m*params.g/params.k;
dl     = 0;
alpha  = 0;
dalpha = 0;
% state vector: x_Stance = [alpha l dalpha dl]
x_sb  = [alpha;l;dalpha;dl];
% energy level
H_bar = H_S(x_sb,params);

x0  = [alpha;l;dalpha;dl];
n   = 4;

% auxiliary variables
xi  = 0; % design parameter

T = 2*pi/sqrt(params.k/params.m);  
t1 = T/2;
t2 = 0;
t3 = T/2;
x1 = x0;
x2 = [params.l_0;0;0;0;0];
x3 = delta_td(x2, params); % x2->x3 in 0 time

u_init = [t3; x3; t2; x2; t1; x1; xi; H_bar];

%% set up continuation
res = @(u) resFun(u,params,ode);

opts             = [];
opts.Grad1       = analyticGradients;
opts.aimOnTarget = false;
opts.idxConPar   = numel(u_init);
opts.FiniteDifferenceStepSize = FiniteDifferenceStepSize;
opts.FunctionTolerance        = FunctionTolerance;
opts.MaxIterations     = 400;
opts.StepTolerance     = 1e-11;
[u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);

HTarget = 2;
h = 5e-2;
d = +1; % not necessary because HTarget is not empty

if u_init(3)==1 && u_init(6)==0
    % apply Lyapunov-Schmidt reduction to branch off into a normal orbit
    tangVec = BifurcationEquation(u_init,res,FiniteDifferenceStepSize);
    idx = find(abs(tangVec(6,:))>FunctionTolerance); % get direction into normal orbit
    u_init = u_init + h*tangVec(:,idx)*sign(tangVec(6,idx));
    [u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);
end

con = numericalContinuation(u_init,jac_init,h,d,HTarget,res,opts);

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
plot(con.u(end,:),con.u(1,:)+con.u(2+n,:))
xlabel('$\bar{H}$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
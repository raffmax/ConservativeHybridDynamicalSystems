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
params.l_0   = 1;   % length of block [units don't matter, could be m]
params.g     = 1;   % gravity [units don't matter, could be m/s^2]
params.r     = 0.01;  % radius of gyration [units of mg/l_0]

% state vector: x = [y phi dy dphi]

% u_init = [t3; x3; t2; x2; t1; x1; xi; H_bar];
u_init = zeros(17,1);

%% set up continuation
res = @(u) resFun(u,params,ode);

opts             = [];
opts.Grad1       = analyticGradients;
opts.aimOnTarget = false;
opts.idxConPar   = numel(u_init);
opts.FiniteDifferenceStepSize = FiniteDifferenceStepSize;
opts.FunctionTolerance        = FunctionTolerance;
opts.MaxIterations     = 10000;
opts.StepTolerance     = 1e-11;
[u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);

if norm(u_init)==0 % check if u_init is the equilibrium
   [~,~,nullVec] = svd(jac_init);
   nullVec = nullVec(:,15:17);
   idx = find(abs(nullVec(16,:))<FunctionTolerance); % remove xi direction
   nullVec = nullVec(:,idx);
   % exploit symmetry: 2*t1 = t2
   coeff_sym = [2*nullVec(1,1)-nullVec(6,1), 2*nullVec(1,2)-nullVec(6,2); nullVec(1,1), nullVec(1,2)]\[0;1];
   tangVec = nullVec*coeff_sym;
   tangVec = tangVec/norm(tangVec);
   u_init = u_init + 0.1*tangVec;
   [u_init,~,~,output,jac_init] = NewtonsMethod(res,u_init,opts);
end

HTarget = 0.6;
h = 2e-2; % fixed step-size
d = +1; % direction of arclength continuation

con = numericalContinuation(u_init,jac_init,h,d,HTarget,res,opts);

%% plots
figure
grid on
hold on
plot(con.u(end,:),con.detSim*1000)
plot(con.u(end,:),con.detAug)
xlabel('$\bar{H}$','Interpreter','latex')
legend('simple Jacobian','augmented Jacobian','Interpreter','latex')

figure
grid on
hold on
plot(con.u(end,:),con.u(1,:)+con.u(6,:)+con.u(11,:))
xlabel('$\bar{H}$','Interpreter','latex')
ylabel('$T$','Interpreter','latex')
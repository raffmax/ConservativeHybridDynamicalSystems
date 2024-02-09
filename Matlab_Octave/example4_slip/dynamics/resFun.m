function [r,dr,Phi_T] = resFun(u,params,ode,varargin)

if nargin < 4
    computePhi_T = false;
else
    computePhi_T = varargin{1};
end

nS = 4;
nF = 5;

% u = [t3; x3; t2; x2; t1; x1; xi; H_bar];
idxOS  = 0; % index offset
t3     = u(idxOS + 1);
idxOS  = idxOS + 1;
x3_bar = u(idxOS + (1:nS));
idxOS  = idxOS + nS;
t2     = u(idxOS + 1);
idxOS  = idxOS + 1;
x2_bar = u(idxOS + (1:nF));
idxOS  = idxOS + nF;
t1     = u(idxOS + 1);
idxOS  = idxOS + 1;
x1_bar = u(idxOS + (1:nS));
idxOS  = idxOS + nS;
xi     = u(idxOS + 1);
idxOS  = idxOS + 1;
H_bar  = u(idxOS + 1);

% simulate phase flows
n = nS;
f1_ode = @(t,x) f_S(x,xi,params);
X0 = [x1_bar;reshape(eye(n),[n^2,1]);zeros(n,1)]; 
if t1==0
    Xt1=X0;
else
    [~, Xt] = ode.solver(f1_ode, [0, t1], X0, ode.settings);
    Xt1 = Xt(end,:)';
end
x1   = Xt1(1:n);
Phi1 = reshape(Xt1(n+(1:n^2)),[n,n]);
psi1 = Xt1(n+n^2+(1:n));
f1   = f1_ode([],x1);

n = nF;
f2_ode = @(t,x) f_F(x,xi,params);
X0 = [x2_bar;reshape(eye(n),[n^2,1]);zeros(n,1)];
if t2==0
    Xt2=X0;
else
    [~, Xt] = ode.solver(f2_ode, [0, t2], X0, ode.settings);
    Xt2 = Xt(end,:)';
end
x2   = Xt2(1:n);
Phi2 = reshape(Xt2(n+(1:n^2)),[n,n]);
psi2 = Xt2(n+n^2+(1:n));
f2   = f2_ode([],x2);

n = nS;
f3_ode = @(t,x) f_S(x,xi,params);
X0 = [x3_bar;reshape(eye(n),[n^2,1]);zeros(n,1)]; 
if t3==0
    Xt3=X0;
else
    [~, Xt] = ode.solver(f3_ode, [0, t3], X0, ode.settings);
    Xt3 = Xt(end,:)';
end
x3   = Xt3(1:n);
Phi3 = reshape(Xt3(n+(1:n^2)),[n,n]);
psi3 = Xt3(n+n^2+(1:n));
f3   = f3_ode([],x3);

% compute reset maps
[x1_P,D12] = delta_lo(x1, params);
[x2_P,D23] = delta_td(x2, params);

% compute events and anchor
[e12, ~, ~, de12] = e_lo(x1, params);
[e23, ~, ~, de23] = e_td(x2, params);
[a, ~, ~, da]     = anchor(x3, params);

% compute first integral
[H1,dH1,~] = H_S(x1_bar,params);

r  = [x3-x1_bar;
      a;
      x2_P-x3_bar;
      e23;
      x1_P-x2_bar;
      e12;
      H1-H_bar];

if computePhi_T
    % compute saltation matrices for monodromy matrix
    f2P = f2_ode([],x2_bar);
    S12 = D12 + (f2P-D12*f1)/(de12*f1);
    
    f3P = f3_ode([],x3_bar);
    S23 = D23 + (f3P-D23*f2)/(de23*f2);
    
    Phi_T = Phi3*S23*Phi2*S12*Phi1;
else
    Phi_T = [];
end

% compute Jacobian
R3  = [f3,Phi3;da*f3,da*Phi3];
R2  = [D23*f2,D23*Phi2;de23*f2,de23*Phi2];
R1  = [D12*f1,D12*Phi1;de12*f1,de12*Phi1];
Xi3 = [psi3;da*psi3];
Xi2 = [D23*psi2;de23*psi2];
Xi1 = [D12*psi1;de12*psi1];

CF  = [zeros(nF,1),-eye(nF);0,zeros(1,nF)];
CS  = [zeros(nS,1),-eye(nS);0,zeros(1,nS)];
OSF = zeros(nS+1,nF+1);   % zero matrix
OFS = zeros(nF+1,nS+1);   % zero matrix
OSS = zeros(nS+1,nS+1);   % zero matrix
%OFF = zeros(nF+1,nF+1);   % zero matrix
OcF  = zeros(nF+1,1); % zero column
OcS  = zeros(nS+1,1); % zero column
OrF  = zeros(1,nF+1); % zero row
OrS  = zeros(1,nS+1); % zero row
h1  = [0,dH1];

dr = [R3  OSF CS  Xi3 OcS;
      CS  R2  OSS Xi2 OcS;
      OFS CF  R1  Xi1 OcF;
      OrS OrF h1  0   -1];
end


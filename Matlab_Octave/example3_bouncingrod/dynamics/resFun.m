function [r,dr,Phi_T] = resFun(u,params,ode,varargin)

if nargin < 4
    computePhi_T = false;
else
    computePhi_T = varargin{1};
end

n = 4;
% u = [t3; x3; t2; x2; t1; x1; xi; H_bar];
idxOS  = 0; % index offset
t3     = u(idxOS + 1);
idxOS  = idxOS + 1;
x3_bar = u(idxOS + (1:n));
idxOS  = idxOS + n;
t2     = u(idxOS + 1);
idxOS  = idxOS + 1;
x2_bar = u(idxOS + (1:n));
idxOS  = idxOS + n;
t1     = u(idxOS + 1);
idxOS  = idxOS + 1;
x1_bar = u(idxOS + (1:n));
idxOS  = idxOS + n;
xi     = u(idxOS + 1);
idxOS  = idxOS + 1;
H_bar  = u(idxOS + 1);

% simulate phase flows
f1_ode = @(t,x) f(x,xi,params);
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

f2_ode = @(t,x) f(x,xi,params);
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

f3_ode = @(t,x) f(x,xi,params);
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
[x1_P,D12] = delta_12(x1, params);
[x2_P,D23] = delta_23(x2, params);

% compute events and anchor
[e12_val, ~, ~, de12] = e_12(x1, params);
[e23_val, ~, ~, de23] = e_23(x2, params);
[a, ~, ~, da]     = anchor(x3, params);

% compute first integral
[H1,dH1,~] = H(x1_bar,params);

r  = [x3-x1_bar;
      a;
      x2_P-x3_bar;
      e23_val;
      x1_P-x2_bar;
      e12_val;
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

C  = [zeros(n,1),-eye(n);0,zeros(1,n)];
Om = zeros(n+1);   % zero matrix
Oc = zeros(n+1,1); % zero column
Or = zeros(1,n+1); % zero row
h1 = [0,dH1];

dr = [R3 Om C  Xi3 Oc;
      C  R2 Om Xi2 Oc;
      Om C  R1 Xi1 Oc;
      Or Or h1  0  -1];
end


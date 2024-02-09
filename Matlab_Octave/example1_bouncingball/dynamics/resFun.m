function [r,dr] = resFun(u,params,ode)

n = 2;

% u = [t2; x2; t1; x1; xi; H_bar];
idxOS  = 0; % index offset
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
f1 = @(t,x) f(x,xi,params);
X0 = [x1_bar;reshape(eye(n),[n^2,1]);zeros(n,1)]; 
if t1==0
    Xt1=X0;
else
    [~, Xt] = ode.solver(f1, [0, t1], X0, ode.settings);
    Xt1 = Xt(end,:)';
end
x1   = Xt1(1:n);
Phi1 = reshape(Xt1(n+(1:n^2)),[n,n]);
psi1 = Xt1(n+n^2+(1:n));
f1   = f1([],x1);

f2 = @(t,x) f(x,xi,params);
X0 = [x2_bar;reshape(eye(n),[n^2,1]);zeros(n,1)]; 
if t2==0
    Xt2=X0;
else
    [~, Xt] = ode.solver(f2, [0, t2], X0, ode.settings);
    Xt2 = Xt(end,:)';
end
x2   = Xt2(1:n);
Phi2 = reshape(Xt2(n+(1:n^2)),[n,n]);
psi2 = Xt2(n+n^2+(1:n));
f2   = f2([],x2);

% compute reset maps
[x1_P,D12] = delta_td(x1, params);

% compute events and anchor
[e12, ~, ~, de12] = e_td(x1, params);
[a, ~, ~, da]     = anchor(x2, params);

% compute first integral
[H1,dH1,~] = H(x1_bar,params);

r  = [x2-x1_bar;
      a;
      x1_P-x2_bar;
      e12;
      H1-H_bar];

% compute Jacobian
R2  = [f2,Phi2;da*f2,da*Phi2];
R1  = [D12*f1,D12*Phi1;de12*f1,de12*Phi1];
Xi2 = [psi2;da*psi2];
Xi1 = [D12*psi1;de12*psi1];
C   = [zeros(n,1),-eye(n);0,zeros(1,n)];
% Om  = zeros(n+1);   % zero matrix
Or  = zeros(1,n+1); % zero row
Oc  = zeros(n+1,1); % zero column
h1  = [0,dH1];

dr = [R2 C  Xi2 Oc;
      C  R1 Xi1 Oc;
      Or h1  0  -1];
end


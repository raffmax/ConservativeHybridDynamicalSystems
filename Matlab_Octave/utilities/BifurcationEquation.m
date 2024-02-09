function [tangVec,Hessian] = BifurcationEquation(u_BP,resFun,FiniteDifferenceStepSize)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[~,dH] = resFun(u_BP);

[F,S,E] = svd(dH);

% dim = sum(diag(S)<FiniteDifferenceStepSize); % dimension of bifurcation
dim = 1;

E1 = E(:,end-dim:end);
E2 = E(:,1:end-dim-1);
F1 = F(:,end-dim+1:end);
F2 = F(:,1:end-dim);

% central differences
ddH = zeros([size(dH),dim+1]);
h = FiniteDifferenceStepSize;
for i = 1:dim+1
    [~,dHdui_P] = resFun(u_BP+h/2*E1(:,i));
    [~,dHdui_N] = resFun(u_BP-h/2*E1(:,i));
    ddH(:,:,i)  = (dHdui_P-dHdui_N)/h;
end

Hessian = zeros(dim+1);
for i = 1:dim+1
    for j = 1:dim+1
        Hessian(i,j) = F1(:,1)'*ddH(:,:,i)*E1(:,j);
    end
end

% g = @(ii,jj) F1'*resFun(u_BP+ii*E1(:,1)+jj*E1(:,2));
% 
% alpha_11 = (g(h,0)-2*g(0,0)+g(-h,0))/(h^2);
% alpha_22 = (g(0,h)-2*g(0,0)+g(0,-h))/(h^2);
% alpha_12 = (g(h,h)+g(-h,-h)-g(h,-h)-g(-h,h))/(4*h^2);
% 
% Hessian_ = [alpha_11,alpha_12;alpha_12,alpha_22];
% 
% [~,jact1] = getJacobianFD(@(u) dhTang(u,resFun,E1(:,1),F1),u_BP,h);
% [~,jact2] = getJacobianFD(@(u) dhTang(u,resFun,E1(:,2),F1),u_BP,h);
% Hessian_2 = [jact1*E1(:,1),jact1*E1(:,2);jact2*E1(:,1),jact2*E1(:,2)];

% degenerate conic section
[V,D] = eig(Hessian);
d = diag(D);
x1 = 1;
y1 = sqrt(-d(1)/d(2));
x2 = 1;
y2 = -y1;

tang1 = [x1;y1];
tang1 = V*tang1/norm(tang1);
tang2 = [x2;y2];
tang2 = V*tang2/norm(tang2);
tangVec = [E1*tang1,E1*tang2];

end

% function dH = dhTang(u,resFun,tang,F1) 
%     [~,jac] = resFun(u);
%     dH = F1'*jac*tang;
% end

function con = numericalContinuation(u_init,jac_init,h,d,HTarget,resFcn,opts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tang_init = getTangent(jac_init);

% Pre-allocate continuation data
con.u      = zeros(numel(u_init),opts.MaxIterations);
con.tang   = zeros(numel(u_init),opts.MaxIterations);
con.detSim = zeros(1,opts.MaxIterations);
con.detAug = zeros(1,opts.MaxIterations);

% Store initial point in continuation data
con.u(:,1)    = u_init;
con.tang(:,1) = tang_init;
con.detSim(1) = det(jac_init(:,1:end-1));
con.detAug(1) = det([jac_init;tang_init']);

% initialize continuation
u    = u_init;
tang = tang_init;

if ~isempty(HTarget)
    dirCon = sign(HTarget-u_init(end));
    % define correct direction of continuation
    if sign(HTarget-u_init(end))*sign(tang_init(end))==-1
        d = -1;
    end
else
    dirCon = 0;
end

endLoop = false;
iter = 1;
while ~endLoop
    % predictor
    u_pred = u+d*h*tang;
    if dirCon*u_pred(end)>dirCon*HTarget
        endLoop          = true;
        opts.aimOnTarget = true;
        u_pred = u_pred-(u_pred(end)-HTarget)*tang/tang(end);
    end
    % corrector
    [u_corr,fval,~,output,jac] = NewtonsMethod(resFcn,u_pred,opts);
    tang_new = getTangent(jac);
    % check for bifurcation
    if tang_new'*tang<0
        display(['Simple Bifurcation between H=',num2str(u(end)),' and H=',num2str(u_pred(end))])
        d = -d;
    elseif tang_new(end)*tang(end)<0
        display(['Turning Point between H=',num2str(u(end)),' and H=',num2str(u_pred(end))])
    end

    % update
    u    = u_corr;
    tang = tang_new;

    if dirCon*u_corr(end)>dirCon*HTarget
        % the corrector found a solution with H>HTarget
        % however, the predictor was H<HTarget
    else
        iter = iter + 1;
        % store continuation data
        con.u(:,iter)    = u;
        con.tang(:,iter) = tang;
        con.detSim(iter) = det(jac(:,1:end-1));
        con.detAug(iter) = det([jac;tang']);
    end

    if endLoop == true
        % iter < opts.MaxIterations
        con.u      = con.u(:,1:iter);
        con.tang   = con.tang(:,1:iter);
        con.detSim = con.detSim(1:iter);
        con.detAug = con.detAug(1:iter);
    end

    % terminate when max iterations
    if iter==opts.MaxIterations
        endLoop = true;
    end
end
end
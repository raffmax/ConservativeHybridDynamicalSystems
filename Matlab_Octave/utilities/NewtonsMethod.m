function[x,fval,exitflag,output,jacobian] = NewtonsMethod(fun,x,options,varargin)
%NewtonsMethod implements Newton's method.
%   NewtonsMethod can be called by the corrector-step of a continuation
%   method or an arbitrary root-finding-problem
%   see also getJacobianFD

if nargin>3
    f        = varargin{1};
    jacobian = varargin{2};
    t        = varargin{3}; % nullspace of jacobian / tangent space
    
    funCounter = 0; % counter for function evaluations
    
    if options.Grad1
        useFD = false; % use finite differences
    else
        useFD    =  true;
        stepSize = options.FiniteDifferenceStepSize;
        % get dimensions of jacobian
        nF    = length(f);
        nX    = length(x);
    end
    
else
   if options.Grad1
       useFD = false;
       [f,jacobian] = fun(x);
       funCounter   = 1; % counter for function evaluations
   else
       useFD        = true;
       stepSize     = options.FiniteDifferenceStepSize;
       [f,jacobian] = getJacobianFD(fun,x,stepSize);
       % get dimensions of jacobian
       nF = length(f);
       nX = length(x);
       funCounter   = 1 + nX*nF; % counter for function evaluations
   end  
   t = getTangent(jacobian);
end

err  = max(abs(f));
iter = 0;
while err > options.FunctionTolerance && iter < options.MaxIterations  
   if err>1
       break
   end
   % newton step
   if options.aimOnTarget
        idx = options.idxConPar;
        delta_x = -[jacobian;[zeros(1,idx-1),1,zeros(1,length(x)-idx)]]\[f;0]; % fix lambda
   else
        delta_x = -[jacobian;t']\[f;zeros(size(t,2))]; % minimze distance to curve if size(t,2)=1 / equivalent to Moore-Penrose Inverse
   end

   if norm(delta_x,1) < options.StepTolerance 
       warning(['The Newton step is smaller than StepTolerance ',num2str(options.StepTolerance),'!'])
       break
   end

   % Newton update
   x = delta_x + x;
   if useFD
       [f,jacobian] = getJacobianFD(fun,x,stepSize);
       funCounter               = funCounter + 1 + nX*nF; 
   else
       [f,jacobian] = fun(x);
       funCounter               = funCounter + 1;
   end
   
   % compute new tangent vector
   t = getTangent(jacobian);

   iter = iter + 1;
   err  = max(abs(f));       
end  
fval   = f;
if err < options.FunctionTolerance
   exitflag = 1;
else
   exitflag = -1;
end

output.t          = getTangent(jacobian);
output.iterations = iter;
output.funcCount  = funCounter;

end

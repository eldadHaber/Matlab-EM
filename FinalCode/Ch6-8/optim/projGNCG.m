function [mc,d] = projGNCG(mc,param)
%[mc,dpre] = projNewtonCG(mc,param)
%  

maxIter     = param.maxIter;     %  Max. no. iterations.
maxCG       = param.maxCGiter;   %  Max cg iters.
cgTol       = param.cgTol;       %  CG stopping tolerance.
stepTol     = param.stepTol;     %  Step norm stopping tol.

low   = param.bounds.low;
high  = param.bounds.high;


regpar = param.regPar;
Wd     = param.Wd;

%---------------------------------------------------------------------------
%  Initialization.
%---------------------------------------------------------------------------
maxStepSize = param.maxStep;

Active = [mc<=low | mc>=high];  % Compute active set

%% evaluate function and derivatives
[sig,dsig] = param.sigFun(mc,param); Jxm = sdiag(dsig);
[d,e]      = param.solveForwardProblem(sig,param);
r          = Wd *(d(:) - param.dobs(:));
[F,dF]     = param.misFun(r,param);
dF         = Jxm'*param.JTmatVec(Wd'*dF,e,sig,param);
    
% compute model objective function
[R,dR,d2R] = param.regFun(mc,param);
    
% objective function
Jc  = F  + regpar*R;
gc  = dF + regpar*dR;

F0 = F; J0 = Jc;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Outer iteration.                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iter = 0;
  outerFlag = 0;
  
  while outerFlag == 0;

    iter = iter + 1;
    fprintf('%3d.0   %3.2e   %3.2e   %3.2e   %3.2e   %3d\n',...
             iter,F/F0,R,regpar,Jc/J0,nnz(Active));
    
    %  Set up Hessian and preconditioner.
    Jx       = @(x) Wd*param.JmatVec(Jxm*x,e,sig,param);
    JTx      = @(x) Jxm*param.JTmatVec(Wd*x,e,sig,param);
    Hs       = @(x) JTx(Jx(x)) + param.regPar*(d2R*(x));
    Precond  = @(x) tril(d2R)\(diag(d2R).*(triu(d2R)\x));
 
    resid = (1-Active) .* gc;
    delm = projPCG(Hs,-resid,Active,Precond,cgTol,maxCG);
    % scale step for the safe size
    if max(abs(delm)) > maxStepSize, delm = delm/max(abs(delm))*maxStepSize; end
    
    % take gradient direction in the active cells
    ga = -gc(mc == low);
    if max(abs(ga)) > max(abs(delm)), ga = ga/max(abs(ga))*max(abs(delm)); end
    delm(mc == low) = ga;
    ga = -gc(mc == high);
    if max(abs(ga)) > max(abs(delm)), ga = ga/max(abs(ga))*max(abs(delm)); end
    delm(mc == high) = ga;
    
    %% Begin line search
    muLS = 1; lsIter = 1;
    while 1,
        mt = mc + muLS*delm;
        mt(mt<low) = low;
        mt(mt>high) = high;
        %% evaluate function and derivatives
        [sig,dsig] = param.sigFun(mt,param); Jxm = sdiag(dsig);
        [d,e]      = param.solveForwardProblem(sig,param);
        r          = Wd *(d(:) - param.dobs(:));
        [F,dF]     = param.misFun(r,param);
    
        % compute model objective function
        [R,dR] = param.regFun(mt,param);
    
        % objective function
        Jt  = F  + regpar*R;
        %%    

        fprintf('%3d.%d   %3.2e   %3.2e   %3.2e   %3.2e\n',iter,lsIter,F/F0,R,regpar,Jt/J0);
        
        if Jt < Jc
            break;
        end
        muLS = muLS/2; lsIter = lsIter+1;
        if lsIter > 6
            disp('LSB '); keyboard
            return
        end
    end
    dF   = Jxm'*param.JTmatVec(Wd'*dF,e,sig,param);
    mnew = mt;
    Jc   = Jt;
    gc = dF + regpar*dR;
    %% End Line search
    
    %% Check for termination 
    stepnorm = norm(mnew-mc);
    mc = mnew;
    Active = [mc<=low | mc>=high];
      
    %  Check stopping criteria for outer iteration.
      
    if iter >= maxIter
      fprintf(' \n\n\n  *** Max iterations emceeded ***\n\n\n\n');
      return
    elseif stepnorm < stepTol * norm(mc(:))
      fprintf('  \n\n\n    *** Step norm tolerance met ***\n');
      return
    end
      
  end % while outer_flag == 0


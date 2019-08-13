function delm = projPCG(Hs,resid,Active,Precond,cgTol,maxCG)
% delm = projPCG(H,resid,Active)
% CG with projection to the active set
%

cgiter = 0;
normResid0 = norm(resid);
delm = zeros(length(resid),1);

while 1
    cgiter = cgiter + 1;
    dc = (1-Active).*(Precond(resid));
    rd = resid'*dc;
   
    %  Compute conjugate direction pc.
    if cgiter == 1,
      pc = dc; 
    else
      betak = rd / rdlast;
      pc = dc + betak * pc;
    end

    %  Form product Hessian*pc.
    Hp = Hs(pc);
    Hp = (1-Active).*Hp;

    %  Update delm and residual.
    alphak = rd / (pc'*Hp);
    delm = delm + alphak*pc;
    resid = resid - alphak*Hp;
    rdlast = rd;
 
    if norm(resid)/normResid0 <= cgTol || cgiter == maxCG 
       break;
    end
end

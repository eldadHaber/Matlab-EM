function[Jv] = getJMatVec(z,e,m,param)
% [v] = getJMatVec(u,m,param)
% 

% magnetic permeability
mu = 4*pi*1e-7 * ones(prod(param.nc),1);
% extract matrices
Curl = param.Curl;
Grad = param.Grad;
Af   = param.Af;
Ae   = param.Ae;
An   = param.An;
Iact = param.Iact;
V    = param.V;
dt   = param.dt;

nt = length(dt); 

% inactive cells
sigma = param.mref + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

e  = reshape(e,[],nt+1);
% Solve linear system
lam = e*0;
Jv  = zeros(size(param.obs,1),nt); 
for i=1:length(dt)

    % The linear system to be solved
    Ke = Curl'*Mmuinv*Curl + 1/dt(i)*Msig;
    
    Gzi = (1/dt(i))*sdiag(e(:,i+1)-(i>1)*e(:,i))*Ae'*V*Iact*z;

    rhs = Gzi + 1/dt(i)*Msig*lam(:,i);
    %lam(:,i+1) = Ke\rhs;
    lam(:,i+1) = solveSystem(Ke,Msig,Mmuinvn,Grad,dt(i),rhs);
    
    % compute Jv
    Jvi     = param.obs*Curl*lam(:,i+1);
    Jv(:,i) = -Jvi;        
end

Jv = Jv(:);

return

% function en = solveSystem(Ke,Msig,Mmuinvn,Grad,dt,rhs)
% 
% % setup preconditioner using Aphi system
% STBa = Grad*Mmuinvn*Grad';
% Aap = [Ke + STBa,         1/dt*Msig*Grad; ...
%        1/dt*Grad'*Msig,   1/dt*Grad'*Msig*Grad];
%     
% Map = @(x) tril(Aap)\(diag(Aap).*(triu(Aap)\x));
% P1  = [speye(size(Ke,2)); Grad'];
% P2  = [speye(size(Ke,2)), Grad];
%     
% MM = @(x) P2*(Map(P1*x));
% 
% [en,~] = pcg(Ke,rhs, 1e-11,1000,MM);

%%

% G = [];
% for i=1:length(dt)
%    if i==1, q=0; else, q=1; end 
%    Gi = 1/dt(i)*sdiag(e(:,i+1)-q*e(:,i))*Ae'*V*Iact;
%    G = [G;Gi];    
% end
% 
% A = Ke;
% B = -1/dt(1)*Msig;          
% Ic = speye(nt);
% Im = spdiags(ones(nt,1),-1,nt,nt);
% AA = kron(Ic,A) + kron(Im,B);
% e1 = zeros(nt,1); e1(1) = 1;
% RR = kron(e1,-1/dt(1)*Msig);
% QQ = kron(Ic,param.obs*Curl);
% e0 = e(:,1);
% 
% dd  = -QQ*(AA\(RR*e0));
% ddf = getTimeDomainData(m,param);
% formoderr1 = norm(dd-ddf(:))
% 
% Msig1 = sdiag(Ae'*(V*(sigma+Iact*z)));
% Ke1   = Curl'*Mmuinv*Curl + 1/dt(i)*Msig1;
% A1    = Ke1;
% B1    = -1/dt(1)*Msig1;          
% AA1   = kron(Ic,A1) + kron(Im,B1);
% RR1   = kron(e1,-1/dt(1)*Msig1);
% 
% [ddf,ee2] = getTimeDomainData(m+z,param);
% e0 = ee2(:,1);
% ee1 = -(AA1\(RR1*e0));
% 
% dd1 = -QQ*(AA1\(RR1*e0));
% formoderr2 = norm(dd1-ddf(:))/norm(dd1(:))
% 
% 
% [norm(dd1-dd),norm(dd1-dd + QQ*(AA\(G*z)))]



%norm(AA*reshape(e(:,2:end),[],1) + RR*e(:,1))
%norm(AA1*reshape(e(:,2:end),[],1) + RR1*e(:,1))
%norm(AA1*reshape(e(:,2:end),[],1) + RR1*e(:,1) - G*z)
%
%
% t = AA\(G*z);
% t1 = lam(:,2:end);
% norm(t(:)-t1(:))
% ee = -AA\(RR*e(:,1));
% norm(ee(:) - reshape(e(:,2:end),[],1)) 
% 
% 0 = d(A(m)u) + A(m)*du + d(q(m)) = A(m)*du + G(u)
% du = -A\G


function[JTv] = getJTmatVec(z,e,m,param)
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

nt = length(dt); ne = size(Curl,2);

% inactive cells
sigma = param.mref + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

e   = reshape(e,[],nt+1);
lam = e*0;

s = Curl'*param.obs'*reshape(z,size(param.obs,1),[]);
JTv = 0;
dt = [dt(:);dt(end)];
for i=length(dt)-1:-1:1
    Ke = Curl'*Mmuinv*Curl + 1/dt(i)*Msig;    
    rhs = s(:,i) + 1/dt(i+1)*Msig*lam(:,i+1);
    lam(:,i) = solveSystem(Ke,Msig,Mmuinvn,Grad,dt(i),rhs);
    Gzi = (1/dt(i))*sdiag(e(:,i+1)-(i>1)*e(:,i))*Ae'*V*Iact;
    JTv = JTv - Gzi'*lam(:,i);
end

end

%%
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
% 
% end

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
w    = param.w;

nw = length(w); ne = size(Curl,2);

% inactive cells
sigma = param.mback + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

e  = reshape(e,ne,[],nw);
ns = size(e,2);
% Solve linear system
Jv  = zeros(size(param.obs,1),ns,nw); 
for i=1:nw
    % The linear system to be solved
    Ke = Curl'*Mmuinv*Curl + 1i*w(i)*Msig;
    
    Gzi = zeros(ne,ns);
    for j=1:ns
        eji      = e(:,j,i);
        Gzi(:,j) = 1i*w(i)*(sdiag(eji)*(Ae'*(V*(Iact*z))));
    end
    %lami = Ke\Gzij;
    lami = solveSystem(Ke,Msig,Mmuinvn,Grad,w(i),Gzi);
    
    % compute Jv
    Jvi       = param.obs*lami;
    Jv(:,:,i) = -Jvi;
end

Jv = Jv(:);

return


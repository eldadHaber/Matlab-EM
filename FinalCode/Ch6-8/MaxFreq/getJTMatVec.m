function[JTv] = getJTMatVec(z,e,m,param)
% [v] = getJTMatVec(u,m,param)
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
z  = reshape(z,[],ns,nw);

% Solve linear system
JTv  = 0; 
for i=1:nw
        
    QTr = param.obs'*z(:,:,i);
    % The linear system to be solved
    Ke   = Curl'*Mmuinv*Curl + 1i*w(i)*Msig;
    lami = solveSystem(Ke',Msig,Mmuinvn,Grad,w(i),QTr);
    
    for j=1:ns
        eji     = e(:,j,i);
        GTlamji = (1i*w(i)*sdiag(eji)*Ae'*V*Iact)'*lami(:,j);
        JTv     = JTv - GTlamji;
    end
end

 JTv = real(JTv);
return


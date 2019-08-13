function[d,e,b] = getFreqDomainData(m,param)
% [d,e,b] = getFreqDomainData(m,param)
% Solve the Backward Euler system 
% Curl*e + i*w*b = 0
% Curl'*(Mmuinv*b) - Msig*e = 0
% 
% By eliminating b and solving for e
% b =  - (1/(i*w))*Curl*e
% (Curl'*Mmuinv*Curl + i*omega*Msig)*en = - i*omega*s
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

% inactive cells
sigma = param.mback + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

w   = param.w;

nw = length(w);
ns = size(param.src,2);

e = zeros(size(Curl,2),ns,nw);
if nargout == 3
    b = zeros(size(Curl,1),ns,nw);
end
d = zeros(size(param.obs,1),ns,nw);

% loop over freq
for i=1:length(w)

    % The linear system to be solved
    Ke = Curl'*Mmuinv*Curl + 1i*w(i)*Msig;
    % solve
    rhs = - 1i*w(i)*param.src;
    %e(:,i+1) = Ke\rhs;
    e(:,:,i) = solveSystem(Ke,Msig,Mmuinvn,Grad,w(i),rhs);
    if nargout == 3
        b(:,:,i) =  -1i*w(i)*Curl*e(:,:,i);
    end
    % compute the data
    d(:,:,i) = param.obs*e(:,:,i);
end

return


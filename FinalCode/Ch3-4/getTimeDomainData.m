function[d,e] = getTimeDomainData(m,param)
% [d,e] = getMisfit(m,param)
% Curl'*Mmuinv*Curl*en + 1/dt*Msig*en = 1/dt * Msig*eo

% magnetic permeability
mu = 4*pi*1e-7 * ones(prod(param.nc),1);
% extract matrices
Curl = param.Curl;  Grad = param.Grad;
Af   = param.Af;      Ae   = param.Ae;
An   = param.An;    Iact = param.Iact;
V    = param.V;
% inactive cells (earth+air)
sigma = param.mref + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

% initial conditions
b0   = param.b0; e0   = Msig\(Curl'*Mmuinv*b0);
dt   = param.dt;

e = zeros(size(Curl,2),length(param.dt)+1);
d = zeros(3,length(dt));
e(:,1) = e0;

% time step
for i=1:length(dt)

    % The linear system to be solved
    Ke = Curl'*Mmuinv*Curl + 1/dt(i)*Msig;
    % solve
    rhs = 1/dt(i)*Msig*e(:,i);
    e(:,i+1) = solveSystem(Ke,Msig,Mmuinvn,Grad,dt(i),rhs);
    % compute the data
    d(:,i) = param.obs*Curl*e(:,i+1);
end

function[d,e,b] = getTimeDomainData(m,param)
% [d,e,b] = getTimeDomainData(m,param)
% Solve the Backward Euler system 
% Curl*e(n) + 1/dt*(b(n) - b(n-1)) = 0
% Curl'*(Mmuinv*b(n))   - Msig*e(n) = 0
% 
% By eliminating b and solving for e
% bn = bo - dt*Curl*en
% Curl'*Mmuinv*(1/dt*bo - Curl*en) - 1/dt*Msig*en = 0
% Curl'*Mmuinv*Curl*en + 1/dt*Msig*en = 1/dt * Msig*eo
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
sigma = param.mref + Iact*m;

% setup mass matrices
Mmuinv  = sdiag(Af'*(V*(1./mu)));
Msig    = sdiag(Ae'*(V*sigma));
Mmuinvn = sdiag(An'*V*(1./mu));

% initial conditions
b0   = param.b0;
e0   = Msig\(Curl'*Mmuinv*b0);
dt   = param.dt;

e = zeros(size(Curl,2),length(param.dt)+1);
b = zeros(size(Curl,1),length(param.dt)+1);
d = zeros(3,length(dt));
e(:,1) = e0; b(:,1) = b0;

Ke = Curl'*Mmuinv*Curl + 1/dt(1)*Msig;
%keyboard
fprintf('Start Chol\n')
C  = chol(Ke);
fprintf('Done Chol\n')

% time step
for i=1:length(dt)

    % The linear system to be solved
    % Ke = Curl'*Mmuinv*Curl + 1/dt(i)*Msig;
    % solve
    rhs = 1/dt(i)*Msig*e(:,i);
    e(:,i+1) = (C\((C')\rhs)); % \rhs;
    %e(:,i+1) = solveSystem(Ke,Msig,Mmuinvn,Grad,dt(i),rhs);
    b(:,i+1) = b(:,i) - dt(i)*Curl*e(:,i+1);
    
    % compute the data
    d(:,i) = param.obs*Curl*e(:,i+1);
    fprintf('time %3.2e\n',sum(dt(1:i)))
end

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
% % The rhs
% 
% en = pcg(Ke,rhs, 1e-11,1000,MM);
% 

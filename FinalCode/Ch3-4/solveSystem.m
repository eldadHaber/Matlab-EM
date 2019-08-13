function en = solveSystem(Ke,Msig,Mmuinvn,Grad,dt,rhs)

% setup preconditioner using Aphi system
STBa = Grad*Mmuinvn*Grad';
Aap = [Ke + STBa,         1/dt*Msig*Grad; ...
       1/dt*Grad'*Msig,   1/dt*Grad'*Msig*Grad];
    
Map = @(x) tril(Aap)\(diag(Aap).*(triu(Aap)\x));
P1  = [speye(size(Ke,2)); Grad'];
P2  = [speye(size(Ke,2)), Grad];
    
MM = @(x) P2*(Map(P1*x));
% The rhs

[en,~] = pcg(Ke,rhs, 1e-9,1000,MM);
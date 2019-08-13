function en = solveSystem(Ke,Msig,Mmuinvn,Grad,dt,rhs)
% Solve the maxwell system using a-phi preconditioner
%
% setup preconditioner using Aphi system
STBa = Grad*Mmuinvn*Grad';
Aap = [Ke + STBa,         1/dt*Msig*Grad; ...
       1/dt*Grad'*Msig,   1/dt*Grad'*Msig*Grad];
 
Map = @(x) tril(Aap)\(diag(Aap).*(triu(Aap)\x));
%v = randn(size(Aap,1),1); v = v/norm(v);   
%lmax = norm(Aap*v)*2; lmin = 0+1e-12;   
%Map = @(x) SolChebyshev(Aap,x,x*0,5,lmax,lmin);
P1  = [speye(size(Ke,2)); Grad'];
P2  = [speye(size(Ke,2)), Grad];


MM = @(x) P2*(Map(P1*x));
% The rhs

en = pcg(Ke,rhs, 1e-11,1000,MM);


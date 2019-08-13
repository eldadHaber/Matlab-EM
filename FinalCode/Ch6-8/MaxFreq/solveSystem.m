function en = solveSystem(Ke,Msig,Mmuinvn,Grad,w,rhs)
% Solve the maxwell system in the freq domain
% using a-phi preconditioner
%

% Ke = CTC + 1i*w*Msig
% Ke*a + STBa*a + 1i*w*Msig*Grad*phi
% 

pcgtol = 1e-7;
% setup preconditioner using Aphi system

STBa = Grad*Mmuinvn*Grad';

if sign(full(imag(sum(diag(Ke))))) > 0
    Aap = [Ke + STBa,         1i*w*Msig*Grad; ...
           1i*w*Grad'*Msig,   1i*w*Grad'*Msig*Grad];
else
    Aap = [Ke + STBa,         -1i*w*Msig*Grad; ...
           -1i*w*Grad'*Msig,  -1i*w*Grad'*Msig*Grad];
end

Map = @(x) tril(Aap)\(sdiag(diag(Aap))*(triu(Aap)\x));

P1  = [speye(size(Ke,2)); Grad'];
P2  = [speye(size(Ke,2)), Grad];


MM = @(x) P2*(Map(P1*x));
%en = Ke\rhs;
en = rhs*0;
for i=1:size(rhs,2)
     en(:,i) = bicgstab(Ke,rhs(:,i), pcgtol,100,MM);
end

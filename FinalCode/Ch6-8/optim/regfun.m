function[R,dR,d2R] = regfun(m,param)
%[R,dR,d2R] = regfun(m,param)
%


G = param.Div' * param.V*param.Iact;

g   = sum(G,2);
ind = find(abs(g)<1e-10);
G = G(ind,:);

flag = 2;
tau = 1e-2;
if flag == 2,

    R  = 0.5*(m-param.mref)'*(G'*(G*(m-param.mref))) + 1e-4*(m(:)-param.mref)'*(m(:)-param.mref);
    dR = G'*(G*(m-param.mref)) + 1e-4*(m(:)-param.mref);
    d2R = G'*G +1e-4*speye(length(m));
elseif flag == 1
    Af = param.Af; Af = Af(:,ind);
    v = diag(param.V);
    R = v'*Af*sqrt((G*(m-param.mref)).^2 + tau);
    dR = G' * sdiag((Af'*v)./sqrt((G*(m-param.mref)).^2 + tau)) * G * (m-param.mref);
    d2R = G' * sdiag((Af'*v)./sqrt((G*(m-param.mref)).^2 + tau)) * G;
end
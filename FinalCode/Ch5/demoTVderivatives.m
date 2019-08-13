t = 1e-8;
n1 = 4; n2 = 5; n3 = 6;
D = getNodalGradientMatrix(n1,n2,n3);
% define function and derivatives
C   = @(u,m)(D'*sdiag(1./sqrt((D*u).^2 + t))*D*u - m);
dCdu = @(u,m)(D'*sdiag(1./sqrt((D*u).^2 + t))*D - ...
                    D'*sdiag((D*u).^2./((D*u).^2 + t).^(3/2))*D);

% now test the derivatives
u = randn(size(D,2),1); m = randn(size(D,2),1);
v = randn(size(u));
f = C(u,m);
A = dCdu(u,m); 
for i=1:10
         h = 10^(-i);
	fp = C(u+h*v,m);
	diff1 = norm(fp-f);
	diff2 = norm(fp-f - h*A*v);
	fprintf('%3.2e    %3.2e    %3.2e\n',h,diff1,diff2)      
end

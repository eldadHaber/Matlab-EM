% define function and derivatives
clear all; close all
n1 = 4; n2 = 5; n3 = 6;
D = getNodalGradientMatrix(n1,n2,n3);
S = getEdgeToCellCenterMatrix(n1,n2,n3);
q = randn((n1+1)*(n2+1)*(n3+1),1);

C       = @(u,m)(D'*sdiag(S'*exp(m))*D*u - q);
dCdm    = @(u,m)(D'*sdiag(D*u)*S'*sdiag(exp(m))); 
dCdu    =  @(u,m)(D'*sdiag(S'*exp(m))*D);

% now test the derivatives
u = randn(size(D,2),1); m = randn(size(S,1),1);
v = randn(size(m));
f = C(u,m);
G = dCdm(u,m); 
for i=1:10
    h = 10^(-i);
	fp = C(u,m+h*v);
	diff1 = norm(fp-f);
	diff2 = norm(fp-f - h*G*v);
	fprintf('%3.2e    %3.2e    %3.2e\n',h,diff1,diff2)      
end

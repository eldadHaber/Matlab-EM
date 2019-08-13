function[Grad] = getNodalGradientMatrix(n1,n2,n3)

ddx = @(n) spdiags(ones(n+1,1)*[-1 1],[0,1],n,n+1);

G1 = kron(speye(n3+1),kron(speye(n2+1),ddx(n1))); 
G2 = kron(speye(n3+1),kron(ddx(n2),speye(n1+1))); 
G3 = kron(ddx(n3),kron(speye(n2+1),speye(n1+1)));
% grad on the nodes 
Grad = [G1; G2; G3];
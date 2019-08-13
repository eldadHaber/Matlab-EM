function[Div] = getFaceDivergenceMatrix(n1,n2,n3)

ddx = @(n) spdiags(ones(n+1,1)*[-1 1],[0,1],n,n+1);

D1 = kron(speye(n3),kron(speye(n2),ddx(n1))); 
D2 = kron(speye(n3),kron(ddx(n2),speye(n1))); 
D3 = kron(ddx(n3),kron(speye(n2),speye(n1)));
% DIV from faces to cell-centers 
Div = [D1 D2 D3];
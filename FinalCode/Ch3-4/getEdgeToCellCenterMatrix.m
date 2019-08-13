function[Aec] = getEdgeToCellCenterMatrix(n1,n2,n3)

av  = @(n) spdiags(ones(n+1,1)*[0.5 0.5],[0,1],n,n+1);

A1 = kron(av(n3),kron(av(n2),speye(n1))); 
A2 = kron(av(n3),kron(speye(n2),av(n1))); 
A3 = kron(speye(n3),kron(av(n2),av(n1)));
% average from edge to cell-centers 
Aec = [A1 A2 A3];

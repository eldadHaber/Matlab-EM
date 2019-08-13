function[Anc] = getNodalToCellCenterMatrix(n1,n2,n3)

av  = @(n) spdiags(ones(n+1,1)*[0.5 0.5],[0,1],n,n+1);

Anc = kron(av(n3),kron(av(n2),av(n1)));

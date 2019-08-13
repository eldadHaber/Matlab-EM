function[Curl] = getEdgeCurlMatrix(n1,n2,n3)

ddx = @(n) spdiags(ones(n+1,1)*[-1 1],[0,1],n,n+1);

nfx = (n1+1)*n2*n3; nfy = n1*(n2+1)*n3; nfz = n1*n2*(n3+1); 
nex = n1*(n2+1)*(n3+1); ney = (n1+1)*n2*(n3+1); nez = (n1+1)*(n2+1)*n3; 


Dyz = kron(ddx(n3),kron(speye(n2),speye(n1+1)));
Dzy = kron(speye(n3),kron(ddx(n2),speye(n1+1)));
 
Dxz = kron(ddx(n3),kron(speye(n2+1),speye(n1)));
Dzx = kron(speye(n3),kron(speye(n2+1),ddx(n1)));

Dxy = kron(speye(n3+1),kron(ddx(n2),speye(n1)));
Dyx = kron(speye(n3+1),kron(speye(n2),ddx(n1)));


% curl on the edges
Curl = [sparse(nfx,nex)          Dyz            -Dzy; ...
        -Dxz             sparse(nfy,ney)         Dzx; ...
        Dxy                      -Dyx     sparse(nfz,nez)];
    
    

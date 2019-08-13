function[Curl,Grad,Div,Af,Ae,An,V] = getLinearOps(h1,h2,h3)
%[Curl,Grad,Div,Af,Ae,An] = getLinearOps(h1,h2,h3)
%

ddx = @(n) spdiags(ones(n+1,1)*[-1 1],[0,1],n,n+1);
av  = @(n) spdiags(ones(n+1,1)*[0.5 0.5],[0,1],n,n+1);

n1 = length(h1); n2 = length(h2); n3 = length(h3);

nfx = (n1+1)*n2*n3; 
nfy = (n2+1)*n1*n3; 
nfz = (n3+1)*n2*n1;

nex = n1*(n2+1)*(n3+1); 
ney = n2*(n1+1)*(n3+1); 
nez = n3*(n2+1)*(n1+1);



%% the Divergence
D1 = kron(speye(n3),kron(speye(n2),ddx(n1))); 
D2 = kron(speye(n3),kron(ddx(n2),speye(n1))); 
D3 = kron(ddx(n3),kron(speye(n2),speye(n1)));
% DIV from faces to cell-centers 
Div = [D1 D2 D3];

%% The nodal grasdiet
G1 = kron(speye(n3+1),kron(speye(n2+1),ddx(n1))); 
G2 = kron(speye(n3+1),kron(ddx(n2),speye(n1+1))); 
G3 = kron(ddx(n3),kron(speye(n2+1),speye(n1+1)));
% grad on the nodes 
Grad = [G1; G2; G3];

%% The Curl from edges to faces
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
    
    
%% Mesh geometry and scaling
V = kron(sdiag(h3),kron(sdiag(h2),sdiag(h1)));

F = sdiag([diag(kron(sdiag(h3),kron(sdiag(h2),speye(n1+1)))); ...
           diag(kron(sdiag(h3),kron(speye(n2+1),sdiag(h1)))); ...
           diag(kron(speye(n3+1),kron(sdiag(h2),sdiag(h1))))]);
                
L = sdiag([diag(kron(speye(n3+1),kron(speye(n2+1),sdiag(h1)))); ...
          diag(kron(speye(n3+1),kron(sdiag(h2),speye(n1+1)))); ...
          diag(kron(sdiag(h3),kron(speye(n2+1),speye(n1+1))))]);
      
Div = V\(Div*F); Grad = L\Grad; Curl = F\(Curl*L); 

%% Face average
A1 = kron(speye(n3),kron(speye(n2),av(n1))); 
A2 = kron(speye(n3),kron(av(n2),speye(n1))); 
A3 = kron(av(n3),kron(speye(n2),speye(n1)));
% average from faces to cell-centers 
Af = [A1 A2 A3];

%% Edge average
A1 = kron(av(n3),kron(av(n2),speye(n1))); 
A2 = kron(av(n3),kron(speye(n2),av(n1))); 
A3 = kron(speye(n3),kron(av(n2),av(n1)));
% average from edge to cell-centers 
Ae = [A1 A2 A3];

%% nodal average
An = kron(av(n3),kron(av(n2),av(n1)));
      
% setup EM problem

% loop size
aloop = 10;
% setup the mesh
npad = 10;   nin  = 4;
padxy = 1.2.^[1:npad];   padz  = 1.2.^[1:npad];
h0    = 100;
h1    = h0*[fliplr(padxy),ones(1,nin+1),padxy]; h1 = h1(:);
h2    = h0*[fliplr(padxy),ones(1,nin+1),padxy]; h2 = h2(:);
h3    = h0*[fliplr(padz),ones(1,nin),padz];         h3 = h3(:);
n1 = length(h1); n2 = length(h2); n3 = length(h3);
nc = [n1, n2, n3];

% loop center
x0 = [sum(h1)/2+h0/2,sum(h2)/2+h0/2,sum(h3)/2+h0/2];

%% Compute differential operators
%[Curl,Grad,Div,Af,Ae,An,V] = getLinearOps(h1,h2,h3);

Curl = getEdgeCurlMatrix(n1,n2,n3);
Grad = getNodalGradientMatrix(n1,n2,n3);
Div  = getFaceDivergenceMatrix(n1,n2,n3);
Af   = getFaceToCellCenterMatrix(n1,n2,n3);
Ae   = getEdgeToCellCenterMatrix(n1,n2,n3);
An   = getNodalToCellCenterMatrix(n1,n2,n3);

[V,F,L] = getMeshGeometry(h1,h2,h3);

Div = V\(Div*F); Grad = L\Grad; Curl = F\(Curl*L);

%% setup active cells (cells that are in the earth)
ind  = reshape(1:n1*n2*n3,n1,n2,n3); ind  = ind(:,:,1:n3/2);
Iact = speye(n1*n2*n3); Iact = Iact(:,ind);

%%
% Compute the initial conditon by computing magetic potential
[Xe1,Ye1,Ze1,Xe2,Ye2,Ze2,Xe3,Ye3,Ze3] = getEdgeGrid(h1,h2,h3);
A = magnetostaticsCurrentLoop([Xe1(:),Ye1(:),Ze1(:)], aloop, x0);  a1 = A(:,1);
A = magnetostaticsCurrentLoop([Xe2(:),Ye2(:),Ze2(:)], aloop, x0);  a2 = A(:,2);
A = magnetostaticsCurrentLoop([Xe3(:),Ye3(:),Ze3(:)], aloop, x0);  a3 = A(:,3);
a  = [a1;a2;a3];

b0 = Curl*a;

%% interpolation matrix for receiver
% Assume measure in the center of the loop
xr = x0;
[Xf1,Yf1,Zf1,Xf2,Yf2,Zf2,Xf3,Yf3,Zf3] = getFaceGrid(h1,h2,h3);

Px = interpmat(unique(Xf1(:)),unique(Yf1(:)),unique(Zf1(:)),xr(1),xr(2),xr(3));
Py = interpmat(unique(Xf2(:)),unique(Yf2(:)),unique(Zf2(:)),xr(1),xr(2),xr(3));
Pz = interpmat(unique(Xf3(:)),unique(Yf3(:)),unique(Zf3(:)),xr(1),xr(2),xr(3));
P  = blkdiag(Px,Py,Pz);

%% lump it all in a structure
param.Curl = Curl;
param.Grad = Grad;
param.Af   = Af;
param.Ae   = Ae;
param.An   = An;
param.V    = V;
param.b0   = b0;
param.Iact = Iact;
param.obs  = P;
param.dt   = ones(20,1)*5e-4;
param.mref = 1e-7;
param.nc   = nc;
param.nf   = prod(nc+[1,0,0]) + prod(nc+[0,1,0]) + prod(nc+[0,0,1]); 

% setup conductivity
m = zeros(n1,n2,n3/2) + 1e-2;  % background
m(fix(n1/2)-3:fix(n1/2)+3,fix(n2/2)-3:fix(n2/2)+3, n3/2-5:n3/2-2) = 1; % block
m = m(:);
% get the data 
[d,e] = getTimeDomainData(m,param);
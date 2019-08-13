% setup Airborne EM problem
clear all
close all

% loop size
aloop = 13;

npad = 10;
nin  = 2;
% setup the mesh
padxy = 1.2.^[1:npad];
padz  = 1.2.^[1:npad];
h0    = 20;
h1    = h0*[fliplr(padxy),ones(1,nin+1),padxy];
h2    = h1;
h3    = h0*[fliplr(padz),ones(1,nin),padz]; 
n1 = length(h1); n2 = length(h2); n3 = length(h3);
h1 = h1(:); h2 = h2(:); h3 = h3(:);

nc = [n1, n2, n3];

% loop center
x0 = [sum(h1)/2+h0/2,sum(h2)/2+h0/2,sum(h3)/2+40];

%% Compute differential operators
[Curl,Grad,Div,Af,Ae,An,V] = getLinearOps(h1,h2,h3);

%% setup active cells (cells that are in the earth)
ind  = reshape(1:n1*n2*n3,n1,n2,n3);
ind  = ind(:,:,1:n3/2);
Iact = speye(n1*n2*n3);
Iact = Iact(:,ind);

%%
% Compute the initial conditon by computing magetic potential
[Xe1,Ye1,Ze1,Xe2,Ye2,Ze2,Xe3,Ye3,Ze3] = getEdgeGrid(h1,h2,h3);
A = magnetostaticsCurrentLoop([Xe1(:),Ye1(:),Ze1(:)], aloop, x0);
a1 = A(:,1);
A = magnetostaticsCurrentLoop([Xe2(:),Ye2(:),Ze2(:)], aloop, x0);
a2 = A(:,2);
A = magnetostaticsCurrentLoop([Xe3(:),Ye3(:),Ze3(:)], aloop, x0);
a3 = A(:,3);
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
param.Div  = Div;
param.Af   = Af;
param.Ae   = Ae;
param.An   = An;
param.V    = V;
param.b0   = b0;
param.Iact = Iact;
param.obs  = P;
t = linspace(0,1e-2,100); %[0,logspace(-7,-2,32)];
param.dt   = diff(t); %ones(32,1)*(1e-2)/32;
param.mref = 1e-7;
param.nc   = nc;
param.nf   = prod(nc+[1,0,0]) + prod(nc+[0,1,0]) + prod(nc+[0,0,1]); 


% setup conductivity
m = 1e-2*(1 + 0*rand(n1,n2,n3/2));
%m(11:16,11:16, 8:10) = 1;
m = m(:);

return
% compute the data, electric and magnetic fields
[d,e,b] = getTimeDomainData(m(:),param);

return
% test derivatives

dm = 1e-4*rand(n1,n2,n3/2);
d1 = getTimeDomainData(m(:)+dm(:),param);
Jdm   = getJMatVec(dm(:),e,m(:),param);
fprintf('%3.2e   %3.2e\n', norm(d1(:)-d(:)),norm(d1(:)-d(:) - Jdm(:)))

% setup the Maxwell system with a random source
clear all; close all
n1 = 4; n2 = 5; n3 = 6;
h1 = 10*rand(4,1);
h2 = 10*rand(5,1);
h3 = 10*rand(6,1);

sigma = rand(n1,n2,n3)*0.1;
mu    = ones(n1,n2,n3)*4*pi*1e-7;

% freq
w = 1e1 * 2*pi;

% linear operators
[V,F,L] = getMeshGeometry(h1,h2,h3);
Afc = getFaceToCellCenterMatrix(n1,n2,n3);
Aec = getEdgeToCellCenterMatrix(n1,n2,n3);
Mu = sdiag(Afc'*(V*(1./mu(:))));
Msig = sdiag(Aec'*(V*(sigma(:))));


CURL = getEdgeCurlMatrix(n1,n2,n3);
DIV  = getFaceDivergenceMatrix(n1,n2,n3);
GRAD = getNodalGradientMatrix(n1,n2,n3);

DIV = V\(DIV*F); GRAD = L\GRAD; CURL = F\(CURL*L);

%% The full Maxwell system and sources
A = [-1i*w*Mu    CURL; CURL'     Msig];

% call for the source (random placeholder)
s = randn(size(CURL,2),1);
b = [zeros(size(CURL,1),1); s];

%% The reduced system
Ae = CURL'* sdiag(1./diag(Mu)) * CURL - 1i*w*Msig;
be = -1i*w*s;

%% Test equivalence
ee = Ae\be;

eh = A\b;

h = eh(1:size(CURL,1)); e = eh(size(CURL,1)+1:end);

fprintf('||e_full - e_red|| =  %3.2e\n',norm(ee-e)/norm(e));


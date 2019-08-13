% Setup the CSEM experiment and data
close all
clear all

addpath ../utils
addpath ../MaxFreq
addpath ../optim
addpath ../EMsources

%% setup the mesh and the domain
h1 = [1.85.^[3:-1:1], ones(1,4),1.85.^[1:1:3]]*50;
h2 = [1.85.^[3:-1:1], ones(1,4),1.85.^[1:1:3]]*50;
h3 = [1.85.^[3:-1:1], ones(1,4),1.85.^[1:1:3]]*25;
h1 = h1(:); h2 = h2(:); h3 = h3(:);

Lx = sum(h1); Ly = sum(h2); Lz = sum(h3);

n1 = length(h1); n2 = length(h2); n3 = length(h3);
nc = n1*n2*n3;

[X,Y,Z] = getCellCenteredGrid(h1,h2,h3);

sigWater = 3;
sigma = zeros(n1,n2,n3)+sigWater;
sigma(:,:,1:n3/2) = 1;
sigma(3:5,3:5,3:4) = 1e-3;



%% setup sources locations
[srcx,srcy] = meshgrid(200:100:1100);
srcx = srcx(:); srcy = srcy(:); srcz = srcx*0+400;

% setup receiver location
[recx,recy] = meshgrid(150:100:1150);
recx = recx(:); recy = recy(:); 
recz = recx*0 + 340;

%% Prepare the matrices for source receiver
ns = length(srcx); nr = length(recx);
% setup source dipoles
x = [0;cumsum(h1)];
y = [0;cumsum(h2)];
z = [0;cumsum(h3)];

nex = n1*(n2+1)*(n3+1);
ney = n2*(n1+1)*(n3+1);
nez = n3*(n2+1)*(n1+1);
ne = nex + ney + nez;
Q = zeros(ne,ns);
for i=1:length(srcx)  
    pt = [srcx(i)-50  srcy(i)  srcz(i); ...
          srcx(i)+50  srcy(i)  srcz(i)];      
    Q(:,i) = getSourceFromLineSeg(pt,x,y,z);
end

% setup receivers (3 components)
P = zeros(ne,3*nr);
cnt = 1;
for i=1:length(recx)  
    pt = [recx(i)-50  recy(i)  recz(i); ...
          recx(i)+50  recy(i)  recz(i)];      
    P(:,cnt) = getSourceFromLineSeg(pt,x,y,z);
    cnt = cnt+1;
end
for i=1:length(recx)  
    pt = [recx(i)  recy(i)-50  recz(i); ...
          recx(i)  recy(i)+50  recz(i)];      
    P(:,cnt) = getSourceFromLineSeg(pt,x,y,z);
    cnt = cnt+1;
end
for i=1:length(recx)  
    pt = [recx(i)  recy(i)  recz(i)-50; ...
          recx(i)  recy(i)  recz(i)+50];      
    P(:,cnt) = getSourceFromLineSeg(pt,x,y,z);
    cnt = cnt+1;
end

% active cells
Iact = speye(nc);
ind = find(sigma~=sigWater); 
Iact = Iact(:,ind);

%% Prepare the param structure
[Curl,Grad,Div,Af,Ae,An,V] = getLinearOps(h1,h2,h3);

param.w = 0.1*2*pi;
param.Curl = Curl;
param.Grad = Grad;
param.Div  = Div;
param.Af   = Af;
param.Ae   = Ae;
param.An   = An;
param.V    = V;

param.Iact  = Iact;
param.nc    = nc;
param.mback = sigWater;

param.obs = P';
param.src = Q;

%%
sigAct = Iact'*sigma(:);

% Compute the forward problem
dobs = getFreqDomainData(sigAct,param);
Wd   = ones(numel(dobs),1)*1/mean(abs(dobs(:)));
%% Now invert

% setup params for inversion
solveForwardProblem = @(m,param)getFreqDomainData(m(:),param);
getJMatVec          = @(z,e,m,param)getJMatVec(z,e,m(:),param);
getJTMatVec         = @(z,e,m,param)getJTMatVec(z,e,m(:),param);
sigFun              = @(m,param)sigfun(m(:),param);
regFun              = @(m,param)regfun(m(:),param);
misFun              = @(r,param)misfun(r(:),param);

param.solveForwardProblem = solveForwardProblem;
param.JmatVec             = getJMatVec;
param.JTmatVec            = getJTMatVec;
param.sigFun              = sigFun;    
param.regFun              = regFun;  
param.misFun              = misFun;

param.dobs      = dobs;
param.Wd        = sdiag(Wd);
%param.epsilond  = 1e-5;
param.regPar    = 1e-16;
param.maxCGiter = 10;
param.cgTol     = 1e-2;
param.maxStep   = 0.5;
param.minUpdate = 1e-4;
param.maxIter   = 20; 
param.stepTol   = 1e-4;
param.mref      = sigAct(:)*0;

param.bounds.low  = log(1e-3);
param.bounds.high = log(1);



% Setup the CSEM experiment and data
close all
clear all

addpath ../utils
addpath ../MaxFreq
addpath ../optim
addpath ../EMsources

%% setup the mesh and the domain
h1 = [1.85.^[3:-1:1], ones(1,24),1.85.^[1:1:3]]*50;
h2 = [1.85.^[3:-1:1], ones(1,24),1.85.^[1:1:3]]*50;
h3 = [1.85.^[3:-1:1], ones(1,24),1.85.^[1:1:3]]*25;
h1 = h1(:); h2 = h2(:); h3 = h3(:);

Lx = sum(h1); Ly = sum(h2); Lz = sum(h3);

n1 = length(h1); n2 = length(h2); n3 = length(h3);
nc = n1*n2*n3;

[X,Y,Z] = getCellCenteredGrid(h1,h2,h3);

% setup salt body
stp = @(t) 0.5*(tanh(100*(t-0.5))+1);
sx = 10^5; sy = 10^5; sz = 5*10^5;
s = stp(exp(-(X-1200).^2/sx - (Y-1200).^2/sy - (Z - 50).^2/sz));
montageArray(s)

% layered earth
m1 = 50*(sin(2*pi*X/Lx) .* sin(2*pi*Y/Ly)) + 600;
m2 = 100*(cos(2*pi*X/Lx) .* cos(2*pi*Y/Ly)) + 400;

sigWater = 3;
sigmaL = zeros(n1,n2,n3); 
sigmaL(Z(:)< m1(:)) = 1;
sigmaL(Z(:)> m1(:)) = sigWater;
sigmaL(Z(:)< m2(:)) = 0.5;

% Reservoir
sx = 10^5; sy = 5*10^5; sz = 2e3;
r = stp(exp(-(X-800).^2/sx - (Y-800).^2/sy - (Z - 500).^2/sz));

% put it all together 
sigma = sigmaL;
sigma(s(:)>0.99) = 5e-2;
sigma(r(:)>0.99) = 5e-4;
sigma(Z(:)> m1(:)) = sigWater;


%% setup sources locations
[srcx,srcy] = meshgrid(500:100:1800);
srcx = srcx(:); srcy = srcy(:); srcz = srcx*0+700;

% setup receiver location
[recx,recy] = meshgrid(350:100:2050);
recx = recx(:); recy = recy(:); 
recz = 50*(sin(2*pi*recx/Lx) .* sin(2*pi*recy/Ly)) + 600;


%% make some pictures if needed
%t1 = [0; cumsum(h1)];
%t2 = [0; cumsum(h2)];
%t3 = [0; cumsum(h3)];
%
% for i=1:34
%   subplot(2,2,1);  
%   imagesc(t1,t3,flipud(squeeze(log10(sigma(:,i,:)))'))
%   axis image
%   subplot(2,2,2);  
%   imagesc(t2,t3,flipud(squeeze(log10(sigma(i,:,:)))'))  
%   axis image
%   subplot(2,2,3.5);  
%   imagesc(t1,t2,squeeze(log10(sigma(:,:,i))))  
%   pause(.5) 
% end

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

%% Prepare the param structure
param.w = 0.1*2*pi;

[Curl,Grad,Div,Af,Ae,An,V] = getLinearOps(h1,h2,h3);
param.Curl = Curl;
param.Grad = Grad;
param.Div  = Div;
param.Af   = Af;
param.Ae   = Ae;
param.An   = An;
param.V    = V;

% active cells
Iact = speye(nc);
ind  = find(sigma~=sigWater); 
Iact = Iact(:,ind);
inda = find(sigma == sigWater);
Ina  = speye(nc);
Ina  = Ina(:,inda);

param.Iact  = Iact;
param.nc    = nc;
param.mback = sigWater*Ina*ones(size(Ina,2),1);

param.obs = P';
param.src = Q;

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

return
[mc,d] = projGNCG(pram.mref-1e-6,param);


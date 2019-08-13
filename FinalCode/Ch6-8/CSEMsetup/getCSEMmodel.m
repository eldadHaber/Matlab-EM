%function[sigma] = getCSEMmodel(n)
%[sigma] = getCSEMmodel(n)
%

close all
n = [128,128,32]*2;

%% setup the mesh and the domain
Lx = 5000; Ly = 5000; Lz = 2500;

h1 = ones(1,n(1))*Lx/n(1);
h2 = ones(1,n(2))*Ly/n(2);
h3 = ones(1,n(3))*Lz/n(3);
h1 = h1(:); h2 = h2(:); h3 = h3(:);

[X,Y,Z] = getCellCenteredGrid(h1,h2,h3);

% setup salt body
stp = @(t) 0.5*(tanh(100*(t-0.5))+1);
sx = 5*10^5; sy = 5*10^5; sz = 1.5*10^6;
salt = stp(exp(-(X-Lx/2).^2/sx - (Y-Ly/2).^2/sy - (Z - 400).^2/sz));

figure(1)
montageArray(salt)

% layered earth
m1 = 50*(sin(2*pi*X/Lx) .* sin(2*pi*Y/Ly)) + 1600;
m2 = 100*(cos(2*pi*X/Lx) .* cos(2*pi*Y/Ly)) + 1400;

sigWater = 3;
sigmaL = zeros(n(1),n(2),n(3)); 
sigmaL(Z(:)< m1(:)) = 1;
sigmaL(Z(:)> m1(:)) = sigWater;
sigmaL(Z(:)< m2(:)) = 0.5;

figure(2)
montageArray(sigmaL)

%return

% Reservoir
sx = 3*10^5; sy = 7*10^5; sz = 2e3;
r = stp(exp(-(X-1900).^2/sx - (Y-1900).^2/sy - (Z - 1300).^2/sz));

% put it all together 
sigma = sigmaL;
sigma(salt(:)>0.99) = 5e-2;
sigma(r(:)>0.99) = 5e-4;
sigma(Z(:)> m1(:)) = sigWater;

montageArray(sigma)



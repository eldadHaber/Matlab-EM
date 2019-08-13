function[Xe1,Ye1,Ze1,Xe2,Ye2,Ze2,Xe3,Ye3,Ze3] = getEdgeGrid(h1,h2,h3)
% [Xe1,Ye1,Ze1,Xe2,Ye2,Ze2,Xe3,Ye3,Ze3] = getEdgeGrid(h1,h2,h3)
% 

n1 = length(h1); n2 = length(h2); n3 = length(h3);
% nodel grids
x = zeros(n1+1,1); for i=1:n1, x(i+1) = x(i) + h1(i); end
y = zeros(n2+1,1); for i=1:n2, y(i+1) = y(i) + h2(i); end
z = zeros(n3+1,1); for i=1:n3, z(i+1) = z(i) + h3(i); end

% cell centered grids
xc = x(1:end-1) + h1/2;
yc = y(1:end-1) + h2/2;
zc = z(1:end-1) + h3/2;

[Xe1,Ye1,Ze1] = ndgrid(xc,y,z);
[Xe2,Ye2,Ze2] = ndgrid(x,yc,z);
[Xe3,Ye3,Ze3] = ndgrid(x,y,zc);




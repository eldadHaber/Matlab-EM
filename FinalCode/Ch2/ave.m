function[Av] = ave(n)
% Av = ave(n)
%
Av = spdiags(ones(n+1,1)*[1/2 1/2],[0,1],n,n+1);
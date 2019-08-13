function[G] = ddx(n)
% [G] = ddx(n)
%
G = spdiags(ones(n+1,1)*[-1 1],[0,1],n,n+1);
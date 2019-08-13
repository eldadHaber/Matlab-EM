function[V] = sdiag(v)
% [V] = sdiag(n)
%
V = spdiags(v(:),0,numel(v),numel(v));
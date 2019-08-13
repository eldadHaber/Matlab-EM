function[mis,dmis,d2mis] = misfun(r,param)
% [mis,dmis,d2mis] = misfun(r,param)
% 

mis   = 0.5*real(r(:)'*r(:));
dmis  = r(:);
d2mis = speye(numel(r));
function[sig,dsig] = sigfun(m,param)
% [sigma] = sigfun(m,param)
%

sig  = exp(m); % + param.bounds.low;
dsig = exp(m);
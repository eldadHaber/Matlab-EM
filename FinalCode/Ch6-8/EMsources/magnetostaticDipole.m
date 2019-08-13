function[A] = magnetostaticDipole(X,p,moment)
%[A] = magnetostaticDipole(X,p)
% a = mu0/(4*pi) * cross(m,r)/abs(r)^3
%

mu0 = 4*pi*1e-7;
n = size(X,1);
R = [X(:,1)-p(1) X(:,2)-p(2) X(:,3)-p(3)];
M = repmat(moment,n,1);

absr = sqrt(R(:,1).^2 + R(:,2).^2 + R(:,3).^2);

MxR = cross(M,R,2);

A = zeros(n,3);
A(:,1) = mu0/(4*pi) * MxR(:,1)./(absr.^3);
A(:,2) = mu0/(4*pi) * MxR(:,2)./(absr.^3);
A(:,3) = mu0/(4*pi) * MxR(:,3)./(absr.^3);

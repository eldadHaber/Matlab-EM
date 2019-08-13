function[J] = getApproxSens(X,xs,xr,w)
%[J] = getApproxSens(X,Y,Z,xs,xr)
%

sig = 3;
mu = 4*pi*1e-7;
k   = sqrt(1i*w*mu*sig);


J = zeros(size(xs,1),numel(X(:,1)));

for j=1:size(xs,1)
    
    R  = sqrt((X(:,1)-xs(j,1)).^2 + (X(:,2)-xs(j,2)).^2 + (X(:,3)-xs(j,3)).^2); 
    S  = sqrt((X(:,1)-xr(j,1)).^2 + (X(:,2)-xr(j,2)).^2 + (X(:,3)-xr(j,3)).^2);
    x  = X(:,1) - xs(j,1);
    xi = X(:,1) - xr(j,1);
    fprintf('%3d\n',j)
    U = 1./(2*pi*sig*R.^3).*(-2 + (1i*k*R+1).*exp(1i*k*R) + 3*x.^2./R.^2) + ...
        1./(2*pi*sig*S.^3).*(-2 + (1i*k*S+1).*exp(1i*k*S) + 3*xi.^2./S.^2);
    J(j,:) = U'; 
    %fprintf('\b\b\b') 
    
end
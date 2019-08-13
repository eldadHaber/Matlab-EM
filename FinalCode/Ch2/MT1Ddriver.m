% setup and test the MT problem in 1D
clear all;
close all;

omega = logspace(1,4,32);
mu    = 4*pi*1e-7;

% setup the basic mesh
sig0 = 1e-2;
skin = sqrt(2./omega/mu/sig0(end)); 


for j = 1:length(omega)

  % setup a deep enough mesh
  L  = 3*skin(j); 
  
  h = diff(linspace(0,2*skin(j),2049));
  h = h(:);
  while 1,
      h = [h;h(end)*1.1];
      if sum(h) > L, break; end
  end
  
  n   = length(h);
  sig = ones(n,1) * sig0;
      
  % mesh     
  z = [0; cumsum(h)];
  % cell-centered mesh
  zc = z(1:end-1) + diff(z)/2;
 
  % the linear system

  G    = ddx(n);
  Av   = ave(n);
  Linv = sdiag(1./h);
  Mmu  = sdiag(h./mu);
  Msig = sdiag(Av'*(sig.*h));
  M    = sdiag(Av'*h);

  % set the matrix system
  % [1i*w         d/dz] [b] - [s1]
  % [1/mu d/dz  -sigma] [e] - [s2]

  A = [1i*omega(j)*speye(n)  Linv*G; ...
        -G'*Linv'*Mmu      -Msig];
    
  % seup boundary conditions b(0) = 1  
  bb = -A(2:end,1);
  A  = A(2:end,2:end);

  % solve
  be = A\bb;
  b = [1;be(1:n-1)];
  e = be(n+1:end);
  
  % extract MT data
  d(j) = e(1);
   
      
end

subplot(2,1,1)
semilogx(omega,1./omega(:)*mu .* abs(d(:)).^2)
title('apparent resistivity')
subplot(2,1,2)
semilogx(omega,angle(d))
title('phase')


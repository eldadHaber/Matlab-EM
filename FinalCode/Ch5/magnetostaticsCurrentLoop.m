function A = magnetostaticsCurrentLoop(x, a, x0)
% A = magnetostaticsCurrentLoop(x, a, x0)
% Compute magnetic vector potential of horizontal circular current loop of
% radius a located at the origin.

n = size(x, 1);
if size(x, 2) ~= 3
  error('x must be an n x 3 array');
end

z = x(:, 3) - x0(3);
y = x(:, 2) - x0(2);
x = x(:, 1) - x0(1);
r = sqrt(x.^2 + y.^2);
m = (4 * a * r) ./ ((a + r).^2 + z.^2);
m(m > 1) = 1; % m might be slightly larger than 1 due to rounding errors
              % but ellipke requires 0 <= m <= 1

[K, E] = ellipke(m);

i = r > 0 & m < 1; % 1/r singular at r = 0 and K(m) singular at m = 1

Aphi = zeros(n, 1);
Ax   = zeros(n, 1);
Ay   = zeros(n, 1);
Az   = zeros(n, 1);

% Common factor is (mu * I) / pi with I = 1 and mu = 4e-7 * pi.
Aphi(i) = 4e-7 ./ sqrt(m(i)) .* sqrt(a ./ r(i)) .* ...
  ((1 - m(i) / 2) .* K(i) - E(i));
Ax(i)   = Aphi(i) .* (-y(i) ./ r(i));
Ay(i)   = Aphi(i) .* ( x(i) ./ r(i));

A = [Ax Ay Az];


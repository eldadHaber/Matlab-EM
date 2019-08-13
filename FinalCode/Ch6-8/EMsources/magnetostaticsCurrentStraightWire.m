function A = magnetostaticsCurrentStraightWire(x, p)
% A = magnetostaticsCurrentStraightWire(x, p)
% Compute magnetic vector potential of a unit line current flowing along
% straight segments given in p.

if size(x, 2) ~= 3
  error('x must be an m x 3 array.');
end
if size(p, 2) ~= 3
  error('p must be an n x 3 array.');
end

m = size(x, 1);
n = size(p, 1);

A = zeros(m, 3);

% End point of segment i is start point of segment i+1. Avoid duplicate
% operations:
b   = p(1, :);
xb  = x - repmat(b, m, 1);
nxb = pythag(xb);

% Sum over all line segments
for i = 1:n-1

  a   = b;
  b   = p(i+1, :);

  xa  = xb;
  xb  = x - repmat(b, m, 1);
  nxa = nxb;
  nxb = pythag(xb);

  t   = b - a;
  nt  = pythag(t);
  if nt == 0
    warning('Polygon segment of zero length detected.'); %#ok<WNTAG>
  end
  t   = t / nt;
  txa = xa * t';
  txb = xb * t';

  jn = nxa > 0 & nxb > 0;
  jb = txb > 0;

  j  = jn &  jb & nxa ~= -txa;
  A(j, :) = A(j, :) + log((nxa(j) + txa(j)) ./ (nxb(j) + txb(j))) * t;
  j  = jn & ~jb & nxa ~=  txa;
  A(j, :) = A(j, :) - log((nxa(j) - txa(j)) ./ (nxb(j) - txb(j))) * t;

end

A = A * 1e-7; % factor mu / 4 / pi = 4e-7 * pi / 4 / pi

function y = pythag(x)
n = size(x, 2);
z = max(abs(x), [], 2);
y = z .* sqrt(sum((x./repmat(z, 1, n)).^2, 2));
y(z == 0) = 0;


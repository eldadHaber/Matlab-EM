function[q,qx,qy,qz] = getSourceFromLineSeg(P,x,y,z)
% [q,qx,qy,qz] = getSourceFromLineSeg(p,x,y,z)
% Generates a current source on the mesh.

% allocate memory for the source
nx = length(x); ny = length(y); nz = length(z);
qx = zeros(nx-1,ny,nz);
qy = zeros(nx,ny-1,nz);
qz = zeros(nx,ny,nz-1);

% discretize each line segment
for i=1:size(P,1)-1
    d = P(i+1,:) - P(i,:);
    dr = find(abs(d)>1e-5);
    if dr == 1,
        qt = interpolateLineSegment(P(i:i+1,:),x,y,z);
        qx = qx + sign(d(dr))*qt;
    elseif dr == 2,
        qt = interpolateLineSegment(P(i:i+1,:),x,y,z);
        qy = qy + sign(d(dr))*qt;
    elseif dr == 3,
        qt = interpolateLineSegment(P(i:i+1,:),x,y,z);
        qz = qz + sign(d(dr))*qt;
    end
end

q = [qx(:);qy(:); qz(:)];
 
%%
% Using interpolation  for different directions


%% Finds the cell for the points
function[cellnum] = findCellNumber(p,x,y,z)
%[cellnum] = findCellNumber(p,x,y,z)
%

[~,im] = min(abs(p(1)-x));
cellnum(1) = im - (p(1) - x(im) < 0);

[~,im] = min(abs(p(2)-y));
cellnum(2) = im - (p(2) - y(im) < 0);

[~,im] = min(abs(p(3)-z));
cellnum(3) = im - (p(3) - z(im) < 0);

%%

function[q] = interpolateLineSegment(p,x,y,z)
% interpolation for different directions

nx = length(x); ny = length(y); nz = length(z);

% determine the direction of the line segment
d  = p(2,:) - p(1,:);
dr = find(abs(d)>1e-5);
if p(end,dr) < p(1,dr), p = flipud(p); end

c1 = findCellNumber(p(1,:),x,y,z);
c2 = findCellNumber(p(2,:),x,y,z);

P(1,:) = p(1,:);
if dr == 1, 
    tt = x;
elseif dr == 2,
    tt = y;
elseif dr == 3,
    tt = z;
end
i = 1;    
while tt(c1(dr)+i) <= tt(c2(dr))
    if dr == 1, 
       P(i+1,:) = [tt(c1(1)+i),p(1,2),p(1,3)];
    elseif dr == 2,
       P(i+1,:) = [p(1,1),tt(c1(2)+i),p(1,3)];
    elseif dr == 3,
       P(i+1,:) = [p(1,1),p(1,2),tt(c1(3)+i)];
    end    
    i = i+1;
end
if tt(c2(dr)) ~= p(2,dr)
  P(i+1,:) = p(2,:);
end

% interpolation weights
if     dr == 1, t1 = y; t2 = z; i1 = 2; i2 = 3;
elseif dr == 2, t1 = x; t2 = z; i1 = 1; i2 = 3;
elseif dr == 3, t1 = x; t2 = y; i1 = 1; i2 = 2;
end

d1(1) = t1(c1(i1)+1) - p(1,i1); 
d1(2) = p(1,i1) - t1(c1(i1));
d2(1) = t2(c1(i2)+1) - p(1,i2); 
d2(2) = p(1,i2) - t2(c1(i2));

a     = sum(d1)*sum(d2);
w1    = (d1(1)*d2(1))/a;
w2    = (d1(2)*d2(1))/a;
w3    = (d1(1)*d2(2))/a;
w4    = (d1(2)*d2(2))/a;

e = zeros(1,3); e(dr) = 1; q = zeros([nx,ny,nz]-e);

for i=1:size(P,1)-1
    L = P(i+1,dr) - P(i,dr);
    if dr == 1
        q(c1(1)+i-1,c1(2),c1(3))     = L*w1;
        q(c1(1)+i-1,c1(2)+1,c1(3))   = L*w2;
        q(c1(1)+i-1,c1(2),c1(3)+1)   = L*w3;
        q(c1(1)+i-1,c1(2)+1,c1(3)+1) = L*w4;
    elseif dr == 2
        q(c1(1),c1(2)+i-1,c1(3))     = L*w1;
        q(c1(1)+1,c1(2)+i-1,c1(3))   = L*w2;
        q(c1(1),c1(2)+i-1,c1(3)+1)   = L*w3;
        q(c1(1)+1,c1(2)+i-1,c1(3)+1) = L*w4;
    elseif dr == 3
        q(c1(1),c1(2),c1(3)+i-1)     = L*w1;
        q(c1(1),c1(2)+1,c1(3)+i-1)   = L*w2;
        q(c1(1)+1,c1(2),c1(3)+i-1)   = L*w3;
        q(c1(1)+1,c1(2)+1,c1(3)+i-1) = L*w4; 
    end
        
end 
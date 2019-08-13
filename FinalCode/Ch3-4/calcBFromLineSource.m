function B = calcBFromLineSource(tx_loc,rx_loc)
%Calculates the analytic h-fields from an arbitrary transmitter
%Inputs: tx_loc -> a nSegments+1 by 3 matrix with the location of
%                   transmittor nodes
%        rx_loc -> nStations by 3 matrix of the locations where the fields
%                   should be calculated


nrx = size(rx_loc,1);
n_segments = size(tx_loc,1) - 1;
mu = 4*pi*1e-7;

B = zeros(nrx,3);

for ii = 1:n_segments

    a = tx_loc(ii,:);
    b = tx_loc(ii+1,:);

    ra = rx_loc - repmat(a,nrx,1);
    rb = rx_loc - repmat(b,nrx,1);

    nra = sqrt(ra(:,1).^2 + ra(:,2).^2 + ra(:,3).^2);
    nrb = sqrt(rb(:,1).^2 + rb(:,2).^2 + rb(:,3).^2);

    ba = b - a;
    nba = norm(ba);

    ra_dot_rb = dot(ra',rb')';

    denom = (nba.^2*nra.^2 - ( ra_dot_rb - nra.^2 ).^2 );
    denom(~denom) = eps;
    t1 = ( ra_dot_rb .* (1./nra + 1./nrb) - ( nra + nrb ) ) ./ denom;

    t2 = [(rb(:,2).*ra(:,3)-rb(:,3).*ra(:,2)), ...
          (rb(:,3).*ra(:,1)-rb(:,1).*ra(:,3)), ...
          (rb(:,1).*ra(:,2)-rb(:,2).*ra(:,1))].*t1(:,ones(1,3));

    B = B + t2;

end
B = B/mu;

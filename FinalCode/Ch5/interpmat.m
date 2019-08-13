function[Q] = interpmat(x,y,z,xr,yr,zr)
%
% This function does local linear interpolation
% computed for each receiver point in turn
%
% 
%
% [Q] = linint(x,y,z,xr,yr,zr)
% Interpolation matrix 
%

 
nx = length(x) ;
ny = length(y) ;
nz = length(z) ;

np = length(xr);

Q = spalloc(np,nx*ny*nz,8*np);

for i = 1:np,

       % fprintf('Point %d\n',i); 

        [dd,im] = min(abs(xr(i)-x));
        
        if  xr(i) - x(im) >= 0,  % Point on the left 
                 ind_x(1) = im;
                 ind_x(2) = im+1;
        elseif  xr(i) - x(im) < 0,  % Point on the right
                 ind_x(1) = im-1; 
                 ind_x(2) = im;
       end;
       dx(1) = xr(i) - x(ind_x(1));
       dx(2) = x(ind_x(2)) - xr(i);

        [dd,im] = min(abs(yr(i) - y)) ; 
       if  yr(i) - y(im) >= 0,     % Point on the left
                 ind_y(1) = im;
                 ind_y(2) = im+1;
       elseif  yr(i) -y(im) < 0,  % Point on the right
                 ind_y(1) = im-1;
                 ind_y(2) = im;
       end;
       dy(1) = yr(i) - y(ind_y(1));
       dy(2) = y(ind_y(2)) - yr(i);

        [dd,im] = min(abs(zr(i) - z));
        if  zr(i) -z(im) >= 0,  % Point on the left
                 ind_z(1) = im;
                 ind_z(2) = im+1;
       elseif  zr(i) -z(im) < 0,  % Point on the right
                 ind_z(1) = im-1;
                 ind_z(2) = im;
       end;
       dz(1) = zr(i) - z(ind_z(1)); 
       dz(2) = z(ind_z(2)) - zr(i);      

       dv = (x(ind_x(2)) - x(ind_x(1))) * (y(ind_y(2)) - y(ind_y(1))) * ...
            (z(ind_z(2)) - z(ind_z(1)));

      Dx =  (x(ind_x(2)) - x(ind_x(1)));
      Dy =  (y(ind_y(2)) - y(ind_y(1)));
      Dz =  (z(ind_z(2)) - z(ind_z(1)));  

      % Get the row in the matrix
      v = zeros(nx, ny,nz);

      v( ind_x(1),  ind_y(1),  ind_z(1)) = (1-dx(1)/Dx)*(1-dy(1)/Dy)*(1-dz(1)/Dz);
      v( ind_x(1),  ind_y(2),  ind_z(1)) = (1-dx(1)/Dx)*(1-dy(2)/Dy)*(1-dz(1)/Dz);
      v( ind_x(2),  ind_y(1),  ind_z(1)) = (1-dx(2)/Dx)*(1-dy(1)/Dy)*(1-dz(1)/Dz);
      v( ind_x(2),  ind_y(2),  ind_z(1)) = (1-dx(2)/Dx)*(1-dy(2)/Dy)*(1-dz(1)/Dz);
      v( ind_x(1),  ind_y(1),  ind_z(2)) = (1-dx(1)/Dx)*(1-dy(1)/Dy)*(1-dz(2)/Dz);
      v( ind_x(1),  ind_y(2),  ind_z(2)) = (1-dx(1)/Dx)*(1-dy(2)/Dy)*(1-dz(2)/Dz);
      v( ind_x(2),  ind_y(1),  ind_z(2)) = (1-dx(2)/Dx)*(1-dy(1)/Dy)*(1-dz(2)/Dz);
      v( ind_x(2),  ind_y(2),  ind_z(2)) = (1-dx(2)/Dx)*(1-dy(2)/Dy)*(1-dz(2)/Dz);

     
      Q(i,:) = (v(:))';

end;

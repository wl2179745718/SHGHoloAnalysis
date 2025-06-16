function [x, y, z, x0_index, y0_index, z0_index, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr)

x = ([0:N(1)-1]-N(1)/2)*dr(1); % generate the x vector

y = ([0:N(2)-1]-N(2)/2)*dr(2); % generate the y vector

z = ([0:N(3)-1]-N(3)/2)*dr(3); % generate the z vector

% spatial frequency vectors
dfx = 1/N(1)/dr(1); % fx spacing
fx = ([0:N(1)-1]-N(1)/2)*dfx; % generate the fx vector
dfy = 1/N(2)/dr(2); % fy spacing
fy = ([0:N(2)-1]-N(2)/2)*dfy; % generate the fy vector

% Object region spatial and spatial frequency coordinates

% make matrix versions of x and y
[X,Y] = meshgrid(x,y);

% make matrix versions of fx and fy
[FX,FY] = meshgrid(fx,fy);

x0_index = find(x==0,1); % where x = 0
y0_index = find(y==0,1);
z0_index = find(z==0,1);

end
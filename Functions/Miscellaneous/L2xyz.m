function [N, x,y,z] = L2xyz(L,delta)
% covert the size of a box to the xyz coordinates in 1D

Lx = L(1); Ly = L(2); Lz = L(3);
dx = delta(1); dy = delta(2); dz = delta(3);

% Nx = Lx/dx; Ny = Ly/dy; Nz = Lz/dz;

% We have to make sure the x,y,z include 0,0,0

%x           = dx* (-Nx/2:1:Nx/2-1);      % 1D axis in x
%y           = dy* (-Ny/2:1:Ny/2-1);      % 1D axis in y
%z           = dz* (-Nz/2:1:Nz/2-1);      % 1D axis in z
x = [ flip( 0:-dx:-Lx/2 ) dx:dx:Lx/2 ];
y = [ flip( 0:-dy:-Ly/2 ) dy:dy:Ly/2 ];
z = [ flip( 0:-dz:-Lz/2 ) dz:dz:Lz/2 ];

Nx = size(x,2);
Ny = size(y,2);
Nz = size(z,2);

N = [Nx, Ny, Nz];

end


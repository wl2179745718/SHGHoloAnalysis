clear variables; close all; clc; addpath(genpath('../Functions'));

% Global Parameters
Units

dx          = 7.08353688324377e-07; dy=dx;% pixel size (x,y) in object space (microns)
dz          = 1.28205128205128e-07;       % pixel size (x,y) in object space (microns)
lambda      =  1030*nm;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
%N           = [N_array(ii)-1, N_array(ii)-1, 2^9];                  % lateral pixel dimension 
L = [0.00145000000000000, 0.00145000000000000, 5.00000000000000e-05];
delta = [dx, dy, dz];
%N = round(L./delta);
deltaf = 1./L;     

% The grid


% x-direction spatial grid
Dx = 1450*um; % length of x window
Nx = 2^11;   % number of x sample points
dx = Dx/(Nx-1); % x sample spacing
x = ([0:Nx-1]-Nx/2)*dx; % generate the x vector

% y-direction spatial grid
Dy = Dx; % length of y window
Ny = Nx;   % number of y sample points
dy = Dy/(Ny-1); % y sample spacing
y = ([0:Ny-1]-Ny/2)*dy; % generate the y vector

% spatial frequency vectors
dfx = 1/Nx/dx; % fx spacing
dfy = 1/Ny/dy; % fy spacing

L = [Dx, Dy, 5.00000000000000e-05];
deltaf = 1./L;
Lf = 1./[dx,dy,dz];
[~,fx,fy,~] = L2xyz(Lf,deltaf);

[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

% The sphere

[X,Y] = meshgrid(x,y);

shape = 'sphere';

%nsphere=1.01;%1.2;
%n=[nsphere,n_imm];
%rad = 2;
%MakeSphereHDF5(rad, n, L, X, Y, z, N, delta,k0);

% The incident field

NA_in = 0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);  
U_in=exp(1i*k*sin(phi)*X);

% The Green's function

%dGk = 1;
%Eps = n_imm^2/lambda^2*0.05;
Gdz  = GOlivier(X, Y, dx, dx, dz, k0);%G_kx_ky(fxx,fyy,n_imm,lambda,dz,dGk,Eps);
%G2dz = GOlivier(X, Y, dx, dx, 2*dz, k0);%G_kx_ky(fxx,fyy,n_imm,lambda,2*dz,dGk,Eps);
%G3dz = GOlivier(X, Y, dx, dx, 3*dz, k0);%G_kx_ky(fxx,fyy,n_imm,lambda,3*dz,dGk,Eps);
%G4dz = GOlivier(X, Y, dx, dx, 4*dz, k0);%G_kx_ky(fxx,fyy,n_imm,lambda,4*dz,dGk,Eps);
Pdz = Propagator(n_imm,lambda,fxx,fyy,dz);
%Gamma = Gammadz(n_imm,lambda,fxx,fyy);



%E_MLB_3d_inc=MultiLayerBornv3(shape,dz,1,U_in,Gdz,Pdz);

%E_in_L=exp(1i*k*cos(phi)*L(3)) * exp(1i*k*sin(phi)* X );

%error = norm(E_in_L - E_MLB_3d_inc,2)/norm(E_MLB_3d_inc,2)

U2 = exp(1i*k*cos(phi)*dz) * U_in;



% make matrix versions of fx and fy
[FX,FY] = meshgrid(fx,fy);

myfft = @(x) fftshift(fft2(ifftshift(x)));
myifft = @(x) ifftshift(ifft2(fftshift(x)));

FR = sqrt(FX.^2+FY.^2);
mysqrt = @(dksq,z) real(sqrt(dksq)) + 1i*sign(z)*imag(sqrt(dksq));
angprop = @(z,lambda) exp(1i*2*pi*z*mysqrt((1/lambda)^2-(FR).^2,z));
eprop = @(e,z,lambda) myifft(myfft(e).*angprop(z,lambda));

UF = myifft(myfft(U_in).*Pdz);

UR = eprop(U_in,dz,lambda);

error = norm(U2-UF,2)/norm(U2,2)

error = norm(UF-UR,2)/norm(UF,2)
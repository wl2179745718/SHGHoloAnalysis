clear variables; close all; clc; addpath(genpath('../Functions'));
% unit = m
ps          = .1;                 % pixel size (x,y,z) in object space (microns)
lambda      = 1;                  % central wavelength (microns)  
n_imm       = 1;                % refractive index of immersion media
nsphere=1.02;%1.2;
n=[nsphere,n_imm];
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [1, 1, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.0;%0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
rad = 2;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.25;

[E_sca_Mie] = Mie_plane(X,Y,z(end),phi,k,rad,c/lambda,n_imm^2,1,nsphere^2,1,40);
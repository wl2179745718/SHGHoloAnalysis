clear variables; close all; clc; addpath(genpath('../../Functions'));

dx          = .1;       % pixel size (x,y) in object space (microns)
dz          = .1;       % pixel size (x,y) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
nsphere=1.01;%1.2;
n=[nsphere,n_imm];
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
%N           = [N_array(ii)-1, N_array(ii)-1, 2^9];                  % lateral pixel dimension 
L = [30, 30, 30];
delta = [dx, dx, dz];
N = round(L./delta);
deltaf = 1./L;
NA_in = 0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

%opt = 'Vol';

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

Eps=n_imm^2/lambda^2*0.05;
dGk = 1;

% sphere
rad = 2;
%RI = MakeSphereInRandMed(rad, n, L, delta);
%V=-(k0)^2*((RI).^2-n_imm^2);
%V_inc = 0*V;
%RI_inc = n_imm * ones(size(RI));

MakeSphereHDF5(rad, n, L, delta);

Gdz = G_kx_ky(fxx,fyy,n_imm,lambda,dz,  dGk,Eps);
G2dz= G_kx_ky(fxx,fyy,n_imm,lambda,2*dz,dGk,Eps);
Pdz = Propagator(n_imm,lambda,fxx,fyy,dz);

U_inp=exp(1i*k*sin(phi)*X);
Incident_HDF5(U_inp,N,Pdz);

Nz = N(3);

U1=U_inp;

V1 = sphere_h5read(1,size(fxx),n_imm,k0);
%V2 = sphere_h5read(2,size(fxx),n_imm);

%Uin2=ifft2(Pdz.*(fft2(U1)));
Uin2 = Incident_h5read(2,size(fxx));

Us2 =(ifft2((fft2(U1.*V1)).*(Gdz))*dz);
U2  =Uin2+Us2;

for i=1:Nz
    Vi = sphere_h5read(i,size(fxx),n_imm,k0);
    
end



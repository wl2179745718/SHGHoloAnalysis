clear variables; close all; clc; addpath(genpath('../Functions'));
% unit = m
ps          = 50;       % pixel size (x,y,z) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
nsphere=1.02;%1.2;
n=[nsphere,n_imm];
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [3, 3, 3];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

rad = 2;

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);

N_o = 40;
error = zeros(1,N_o-1);

E_sca_Mie_1 = Mie_plane(X,Y,z(end),phi,k,rad,c/lambda,n_imm^2,1,nsphere^2,1,1);
for ii = 2:N_o

E_sca_Mie = Mie_plane(X,Y,z(end),phi,k,rad,c/lambda,n_imm^2,1,nsphere^2,1,ii);
error(ii-1) = norm(E_sca_Mie-E_sca_Mie_1)/norm(E_sca_Mie_1);

end

plot(error)


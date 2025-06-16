clear variables; close all; clc; addpath(genpath('../Functions'));

dx          = .1;       % pixel size (x,y) in object space (microns)
dz          = .1;       % pixel size (x,y) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
k0=(2*pi)/lambda;
k=k0*n_imm;

Lmax = 300;

L = [Lmax, Lmax, 100];
delta = [dx, dx, dz];
%N = round(L./delta);
%deltaf = 1./L;

[N, x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);  

Eps=0;%n_imm^2/lambda^2*0.05;
dGk = 1; % 1 is no clip, 0 is full clip to the singularity

Gdz = fftshift(G_kx_ky(fxx,fyy,n_imm,lambda,dz,  dGk,Eps));
%Gdz = G_kx_ky(fxx,fyy,n_imm,lambda,dz,  dGk,Eps);

%G = -exp(1i*k0*sqrt(dz^2+ X.^2 + Y.^2))./sqrt(dz^2+ X.^2 + Y.^2);

%Gk = fftshift(ifftshift(fft2(G)))*dx^2/pi/4;
%Gk = myfft(G)*dx^2/pi/4;

%fx = fftshift(fx);

Gk = GOlivier(X, Y, dx, dx, dz, k0);

figure
subplot(2,2,1)
imagesc(fx,fx,abs(Gdz));
set(gca,'YDir','normal')
title('|Gdz|')
axis equal

subplot(2,2,2)
imagesc(fx,fx,abs(Gk));
set(gca,'YDir','normal')
title('|Gk|')
axis equal

subplot(2,2,3)
hold on
plot(abs(Gdz(round(end/2),:)))
plot(abs(Gk(round(end/2),:)))
legend('Gdz','Gk')
subplot(2,2,4)
hold on
plot(angle(Gdz(round(end/2),:)))
plot(angle(Gk(round(end/2),:)))
legend('Gdz','Gk')

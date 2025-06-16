clear variables; close all; clc; addpath(genpath('../../Functions'));

dx          = .1;       % pixel size (x,y) in object space (microns)
dz          = .1;       % pixel size (x,y) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
k0=(2*pi)/lambda;
k=k0*n_imm;

L = [30, 30, 30];
delta = [dx, dx, dz];
N = round(L./delta);
deltaf = 1./L;

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);  

Eps=0;%n_imm^2/lambda^2*0.05;
dGk = 1;

Gdz = fftshift(G_kx_ky(fxx,fyy,n_imm,lambda,dz,  dGk,Eps));

G = exp(1i*k0*sqrt(dz^2+ X.^2 + Y.^2))./sqrt(dz^2+ X.^2 + Y.^2);

Gk = fftshift(fft2(G))/L(1)^2;

fx = fftshift(fx);

figure
subplot(2,2,1)
imagesc(fx,fx,abs(Gdz));
set(gca,'YDir','normal')
title('|Gdz|')
subplot(2,2,2)
imagesc(fx,fx,abs(Gk));
set(gca,'YDir','normal')
title('|Gk|')
subplot(2,2,3)
hold on
plot(fx,abs(Gdz(round(end/2),:)))
plot(fx,abs(Gk(round(end/2),:)))
legend('Gdz','Gk')
subplot(2,2,4)
hold on
plot(angle(Gdz(round(end/2),:)))
plot(angle(Gk(round(end/2),:)))
legend('Gdz','Gk')


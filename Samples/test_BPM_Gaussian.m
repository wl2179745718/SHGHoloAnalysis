clear variables; close all; clc; addpath(genpath('../Functions'));

% length unit is lambda
k = 2*pi;
w0 = 5;
Z  = -100;
zR = k*w0^2/2;

delta = [0.05,0.05,0.05];
L = [50,50,50];
[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

rho = sqrt(X.^2+Y.^2);
E_z = 1/(1+1i*Z/zR)*exp(1i*k*Z)*exp(-rho.^2./(w0^2*(1+1i*Z/zR)));
E_0 = BPM(E_z, 0-Z, inf, k, 2*pi*fxx, 2*pi*fyy);
E_w = exp(-rho.^2./(w0^2));

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,3,1)
imagesc(x,y,abs(E_z))
axis square
%clim([0 cmax]);
colorbar
title('|E_{z}|')

subplot(2,3,4)
imagesc(x,y,angle(E_z))
axis square
colorbar
title('\angle E_{z}')

subplot(2,3,2)
imagesc(x,y,abs(E_0))
axis square
%clim([0 cmax]);
colorbar
title('|E_{0}|')

subplot(2,3,5)
imagesc(x,y,angle(E_0))
axis square
colorbar
title('\angle E_{0}')

subplot(2,3,3)
imagesc(x,y,abs(E_w))
axis square
%clim([0 cmax]);
colorbar
title('|E_{w}|')

subplot(2,3,6)
imagesc(x,y,angle(E_w))
axis square
colorbar
title('\angle E_{w}')

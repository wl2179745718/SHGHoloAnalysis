clear variables; close all; clc; addpath(genpath('../Functions'));
%addpath '/Users/langwang/Documents/GitHub/MultiSlabRytovCode'

ps          = .05;       % pixel size (x,y,z) in object space (microns)
N           = [2^10-1, 2^10-1, 2^9];                  % lateral pixel dimension 
delta = [ps, ps, ps];
L = delta.*N;
[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy
lambda = 1;
n_imm = 1;
k0=(2*pi)/lambda;
k=k0*n_imm;

l=.5;                 % Correlation length (um)
sigma=.01;            % standard deviation of RI
nu = 1;
RImodel='MW';

[RI sdev avg]=RIGenerator3D(RImodel,l,sigma,nu,x,y,z,'single');

RI=RI-avg;
RI=n_imm+RI*(sigma/sdev); % Rescaling std
RI = n_imm+0*RI;
V=-(k0)^2*((RI).^2-n_imm^2);

mu_sMWa=@(l,nu,sig,k) 2*sqrt(pi)*(k^2).*l.*(sig^2).*(gamma(nu+.5)./abs(gamma(nu))).*(1-(1+4.*(k.*l).^2).^(-nu-.5));
mu_sMW=mu_sMWa(l,nu,sigma,k);
L1=0.2*(1/mu_sMW);
NA=n_imm*sin(x(end)*pi*n_imm*.61/2/L1);              % Set NA so that the size of the beam at the entrance of the material is small enough such that it will not be cut off
W0=.61*lambda/asin(NA/n_imm);                      % Beam waist at center wavelength
z0v=pi*n_imm*W0^2/lambda;
Rv=@(z0,z) z.*(1+(z0./z).^2);
wv=@(z0,w0,z) sqrt(w0^2*(1+(z/z0).^2));
Ever=@(z0,w0,k,z) (w0./(wv(z0,w0,z))).*exp(-(sqrt(X.^2+Y.^2).^2)./(wv(z0,w0,z).^2)).*exp(-1i.*(k.*z-atan(z./z0))).*exp(-1i.*((k*sqrt(X.^2+Y.^2).^2)./(2*Rv(z0,z))));
U_in=Norm(Ever(z0v,W0,k,L1));

w0 = 0.6*lambda;
ZR = pi*n_imm/lambda*w0^2;
ZG = -0.8*L(3);
NA_G = n_imm*w0/ZR;
U_in = 1/4/pi/(1i*ZG+ZR)*exp(1i*2*pi*n_imm/lambda*ZG)*exp(-pi*n_imm/lambda/(1i*ZG+ZR)*(X.^2+Y.^2));
Zmed = ZG+0.5*L(3);
U_med= 1/4/pi/(1i*Zmed+ZR)*exp(1i*2*pi*n_imm/lambda*Zmed)*exp(-pi*n_imm/lambda/(1i*Zmed+ZR)*(X.^2+Y.^2));


opt = 'Vol';
ord = 1;
dGk = 1;%-1e-6;
Eps = n_imm^2/lambda^2*0.05;%0%n_imm^2/lambda^2*1e-9;

E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_in,ord,Eps,dGk,opt);
E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_in,Eps,dGk,opt);%.017
E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_in,Eps,dGk,opt);%.017
E_MSL_3d=MultiSlice(U_in,RI,lambda,ps,fxx,fyy,n_imm,opt);

V_inc = 0*V;
E_MLB_3d0=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_in,Eps,dGk,opt);%.017

figure
colormap(gray(256));
imagesc(z,x,squeeze(RI(:,round(end/2),:)))
colorbar

figure
subplot(1,2,1)
imagesc(x,y,abs(U_in))
colorbar
subplot(1,2,2)
imagesc(x,y,angle(U_in))
colorbar

figure
set(gcf, 'Position', get(0, 'Screensize'));
%sgtitle('\epsilon = 0.05 n_{imm}^2/\lambda^2')
colormap(gray(256));
subplot(2,5,1)
imagesc(z,x,abs(squeeze(E_MLB_3d(:,round(end/2),:))));
pbaspect([1 round(L(1)/L(3)) 1])
title('MLB')
subplot(2,5,2)
imagesc(z,x,abs(squeeze(E_MLR_3d(:,round(end/2),:))));
pbaspect([1 round(L(1)/L(3)) 1])
title('MLR')
subplot(2,5,3)
imagesc(z,x,abs(squeeze(E_MSR_3d(:,round(end/2),:))));
pbaspect([1 round(L(1)/L(3)) 1])
title('MSR')
subplot(2,5,4)
imagesc(z,x,abs(squeeze(E_MSL_3d(:,round(end/2),:))));
pbaspect([1 round(L(1)/L(3)) 1])
title('MSL')
subplot(2,5,5)
imagesc(z,x,abs(squeeze(E_MLB_3d0(:,round(end/2),:))));
pbaspect([1 round(L(1)/L(3)) 1])
title('Homogeneous')

subplot(2,5,6)
imagesc(x,y,abs(squeeze(E_MLB_3d(:,:,round(end/2)))));
colorbar
pbaspect([1 1 1])
subplot(2,5,7)
imagesc(x,y,abs(squeeze(E_MLR_3d(:,:,round(end/2)))));
colorbar
pbaspect([1 1 1])
subplot(2,5,8)
imagesc(x,y,abs(squeeze(E_MSR_3d(:,:,round(end/2)))));
colorbar
pbaspect([1 1 1])
subplot(2,5,9)
imagesc(x,y,abs(squeeze(E_MSL_3d(:,:,round(end/2)))));
colorbar
pbaspect([1 1 1])
subplot(2,5,10)
imagesc(x,y,abs(U_med));
colorbar
pbaspect([1 1 1])

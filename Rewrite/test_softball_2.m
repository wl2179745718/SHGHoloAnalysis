clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 10*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [400*um, 400*um, 20*um];
dr  = [ 708*nm,  708*nm,20*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
thetain = 15/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
%shape = 'sphere';
%nsphere=1.02;%1.2;
%rad = 10*um;
%MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

shape = 'softball';
n_center=1.002;%1.2;
sigma = 2.5*um;
MakeSoftballHDF5(sigma, [n_center,n_imm], Box, X, Y, z,N,dr, k0);

% small box for MOM
box = [4*sigma, 4*sigma, 4*sigma];

% Correct ouput wave
uout = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

%uouttst = Propagator(uin,lambda,FX,FY,Box(3));
%PZ = Propagator(lambda,FX,FY,z(end)-z(1));
%uouttst = myifft(myfft(uin).*PZ);

%error = norm(uout-uouttst,2)/norm(uout,2)

% initialize
%udL = uin;
Pdz  = Propagator(lambda,FX,FY,dr(3));
P2dz = Propagator(lambda,FX,FY,2*dr(3));
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);
G2dz = GOlivier(X, Y, dr(1), dr(2), 2*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,2*dz,dGk,Eps);
G3dz = GOlivier(X, Y, dr(1), dr(2), 3*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,3*dz,dGk,Eps);
Gamma = Gammadz(lambda,FX,FY);
%for ii=2:N(3)

    % propagate slice
%    udL = myifft(myfft(udL).*Pdz);

%end

% slice error
%error = norm(uout-udL,2)/norm(uout,2)

U_MLB=MLB(shape,uin,Gdz,Pdz);
U_MLR=MLR(shape,uin,Gdz,Pdz);
orderMSR=2;
U_MSR=MSR(shape,uin,orderMSR,Gdz,Pdz,Gamma);
U_MLB2 = MLB2order(shape,uin,Gdz,Pdz,P2dz);
U_MLB4 = MLB4order(shape,uin,Gdz,G2dz,G3dz,Pdz);

U_MOM = MOM(shape, X, Y, box, z, fxin, thetain, lambda, k0);

norm(U_MLB-U_MLR,2)/norm(U_MLB,2)
norm(U_MLB-U_MSR,2)/norm(U_MLB,2)
norm(U_MLR-U_MSR,2)/norm(U_MLR,2)
norm(U_MLB-U_MLB2,2)/norm(U_MLB,2)
norm(U_MLB-U_MLB4,2)/norm(U_MLB,2)

figure_window = [-15*sigma 15*sigma];

figure
subplot(3,4,1)
imagesc(x,y,abs(U_MLB-uout))
axis square
colorbar
title('|U_{MLB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,2)
imagesc(x,y,angle(U_MLB-uout))
axis square
colorbar
title('\angle U_{MLB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,5)
imagesc(x,y,abs(U_MLB2-uout))
axis square
colorbar
title('|U_{MLB2}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,6)
imagesc(x,y,angle(U_MLB2-uout))
axis square
colorbar
title('\angle U_{MLB2}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,9)
imagesc(x,y,abs(U_MLB4))
axis square
colorbar
title('|U_{MLB4}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,10)
imagesc(x,y,angle(U_MLB4))
axis square
colorbar
title('\angle U_{MLB4}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,3)
imagesc(x,y,abs(U_MLR-uout))
axis square
colorbar
title('|U_{MLR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,4)
imagesc(x,y,angle(U_MLR-uout))
axis square
colorbar
title('\angle U_{MLR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,7)
imagesc(x,y,abs(U_MSR-uout))
axis square
colorbar
title('|U_{MSR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,8)
imagesc(x,y,angle(U_MSR-uout))
axis square
colorbar
title('\angle U_{MSR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,11)
imagesc(x,y,abs(U_MOM))
axis square
colorbar
title('|U_{MOM}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,12)
imagesc(x,y,angle(U_MOM))
axis square
colorbar
title('\angle U_{MOM}')
xlim(figure_window)
ylim(figure_window)

set(gcf, 'Position', get(0, 'Screensize'));

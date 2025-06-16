clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 1*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [40*um, 40*um, 20*um];
dr  = [60*nm, 60*nm, 80*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
thetain = 15/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
shape = 'sphere';
nsphere=1.002;%1.2;
rad = 2*um;
MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

%shape = 'softball';
%n_center=1.002;%1.2;
%sigma = 2*um;
%MakeSoftballHDF5(sigma, [n_center,n_imm], Box, X, Y, z,N,dr, k0);

% small box for MOM
box = [2.2*rad, 2.2*rad, 2.2*rad];

% Correct ouput wave
uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

[X1,Z1] = meshgrid(x,z);
uin_xz = exp(1i*2*pi*fxin*X1).*exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

[Y1,Z1] = meshgrid(y,z);
uin_yz = exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

x0_index = find(x==0,1);
y0_index = find(y==0,1);

% initialize
%udL = uin;
Pdz  = Propagator(lambda,FX,FY,dr(3));
Gamma = Gammadz(lambda,FX,FY);
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);
Diagnose_Green(Gdz, fx, fy, lambda);
[U_MLB_xy, U_MLB_xz, U_MLB_yz]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);


figure_window = [-3*rad/um 3*rad/um];

figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(U_MLB_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(U_MLB_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(U_MLB_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(U_MLB_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(U_MLB_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(U_MLB_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

sgtitle('MLB')
set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',30)


[U_MLR_xy, U_MLR_xz, U_MLR_yz]=MLR(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);


figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(U_MLR_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(U_MLR_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(U_MLR_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(U_MLR_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(U_MLR_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(U_MLR_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

sgtitle('MLR')
set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',30)


order=2;
[U_MSR_xy, U_MSR_xz, U_MSR_yz]=MSR(shape,uin,order,Gdz,Pdz, Gamma, uin_xy, uin_xz, uin_yz, x0_index, y0_index);


figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(U_MSR_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(U_MSR_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(U_MSR_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(U_MSR_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(U_MSR_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(U_MSR_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

sgtitle('MSR')
set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',30)


U_MOM = MOM(shape, X, Y, box, z, fxin, thetain, lambda, k0);


dr = [dr(1), dr(2), dr(3)/2];
[Box,N] = Box_Regularization(Box,dr);
[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
%shape = 'sphere';
%nsphere=1.02;%1.2;
%rad = 10*um;
%MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

[X1,Z1] = meshgrid(x,z);
uin_xz = exp(1i*2*pi*fxin*X1).*exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

[Y1,Z1] = meshgrid(y,z);
uin_yz = exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

x0_index = find(x==0,1);
y0_index = find(y==0,1);

Pdz  = Propagator(lambda,FX,FY,dr(3));
P2dz = Propagator(lambda,FX,FY,2*dr(3));
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);

[U_MLB2_xy, U_MLB2_xz, U_MLB2_yz]=MLB2order(shape,uin,Gdz,Pdz,P2dz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);


figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(U_MLB2_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(U_MLB2_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(U_MLB2_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(U_MLB2_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(U_MLB2_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(U_MLB2_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

sgtitle('MLB 2nd order')
set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',30)


dr = [dr(1), dr(2), dr(3)/2];
[Box,N] = Box_Regularization(Box,dr);
[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
%shape = 'sphere';
%nsphere=1.02;%1.2;
%rad = 10*um;
%MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

Pdz  = Propagator(lambda,FX,FY,dr(3));
Gamma = Gammadz(lambda,FX,FY);
P2dz = Propagator(lambda,FX,FY,2*dr(3));
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);
G2dz = GOlivier(X, Y, dr(1), dr(2), 2*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,2*dz,dGk,Eps);
G3dz = GOlivier(X, Y, dr(1), dr(2), 3*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,3*dz,dGk,Eps);

x0_index = find(x==0,1);
y0_index = find(y==0,1);

[U_MLB4_xy, U_MLB4_xz, U_MLB4_yz] = MLB4order(shape,uin,Gdz,G2dz,G3dz,Pdz, x0_index, y0_index);


dr = [dr(1), dr(2), dr(3)*4];
[Box,N] = Box_Regularization(Box,dr);
[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(U_MLB4_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(U_MLB4_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(U_MLB4_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(U_MLB4_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(U_MLB4_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(U_MLB4_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

sgtitle('MLB 4th order')
set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',30)



figure_window = [-5*rad 5*rad];

figure
subplot(3,4,1)
imagesc(x,y,abs(U_MLB_xy))
axis square
colorbar
title('|U_{MLB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,2)
imagesc(x,y,angle(U_MLB_xy))
axis square
colorbar
title('\angle U_{MLB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,3)
imagesc(x,y,abs(U_MLR_xy))
axis square
colorbar
title('|U_{MLR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,4)
imagesc(x,y,angle(U_MLR_xy))
axis square
colorbar
title('\angle U_{MLR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,7)
imagesc(x,y,abs(U_MSR_xy))
axis square
colorbar
title('|U_{MSR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,8)
imagesc(x,y,angle(U_MSR_xy))
axis square
colorbar
title('\angle U_{MSR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,5)
imagesc(x,y,abs(U_MLB2_xy))
axis square
colorbar
title('|U_{MLB2}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,6)
imagesc(x,y,angle(U_MLB2_xy))
axis square
colorbar
title('\angle U_{MLB2}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,9)
imagesc(x,y,abs(U_MLB4_xy))
axis square
colorbar
title('|U_{MLB4}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,10)
imagesc(x,y,angle(U_MLB4_xy))
axis square
colorbar
title('\angle U_{MLB4}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,11)
imagesc(x,y,abs(U_MOM))
axis square
colorbar
title('|U_{1stB}|')
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



figure
subplot(3,4,1)
imagesc(x,y,abs(U_MLB_xy-U_MOM))
axis square
colorbar
title('|U_{MLB}-U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,2)
imagesc(x,y,angle(U_MLB_xy-U_MOM))
axis square
colorbar
title('\angle U_{MLB}-U_{1stB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,3)
imagesc(x,y,abs(U_MLR_xy-U_MOM))
axis square
colorbar
title('|U_{MLR}-U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,4)
imagesc(x,y,angle(U_MLR_xy-U_MOM))
axis square
colorbar
title('\angle U_{MLR}-U_{1stB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,7)
imagesc(x,y,abs(U_MSR_xy-U_MOM))
axis square
colorbar
title('|U_{MSR}-U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,8)
imagesc(x,y,angle(U_MSR_xy-U_MOM))
axis square
colorbar
title('\angle U_{MSR}-U_{1stB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,5)
imagesc(x,y,abs(U_MLB2_xy-U_MOM))
axis square
colorbar
title('|U_{MLB2}-U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,6)
imagesc(x,y,angle(U_MLB2_xy-U_MOM))
axis square
colorbar
title('\angle U_{MLB2}-U_{1stB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,9)
imagesc(x,y,abs(U_MLB4_xy-U_MOM))
axis square
colorbar
title('|U_{MLB4}-U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,10)
imagesc(x,y,angle(U_MLB4_xy-U_MOM))
axis square
colorbar
title('\angle U_{MLB4}-U_{1stB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,11)
imagesc(x,y,abs(U_MOM))
axis square
colorbar
title('|U_{1stB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,12)
imagesc(x,y,angle(U_MOM))
axis square
colorbar
title('\angle U_{1stB}')
xlim(figure_window)
ylim(figure_window)

set(gcf, 'Position', get(0, 'Screensize'));

lmax=20;
angle_correction = exp(1i*2*pi*(0-z(1)+dr(3))*cos(thetain)/lambda);
U_Mie = angle_correction.*ScalarMieField(X,Y,z,thetain,lmax,n_imm,nsphere,k0,rad);


figure
subplot(3,4,1)
imagesc(x,y,abs(U_MLB_xy))
axis square
colorbar
title('|U_{MLB}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,2)
imagesc(x,y,angle(U_MLB_xy))
axis square
colorbar
title('\angle U_{MLB}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,3)
imagesc(x,y,abs(U_MLR_xy))
axis square
colorbar
title('|U_{MLR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,4)
imagesc(x,y,angle(U_MLR_xy))
axis square
colorbar
title('\angle U_{MLR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,7)
imagesc(x,y,abs(U_MSR_xy))
axis square
colorbar
title('|U_{MSR}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,8)
imagesc(x,y,angle(U_MSR_xy))
axis square
colorbar
title('\angle U_{MSR}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,5)
imagesc(x,y,abs(U_MLB2_xy))
axis square
colorbar
title('|U_{MLB2}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,6)
imagesc(x,y,angle(U_MLB2_xy))
axis square
colorbar
title('\angle U_{MLB2}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,9)
imagesc(x,y,abs(U_MLB4_xy))
axis square
colorbar
title('|U_{MLB4}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,10)
imagesc(x,y,angle(U_MLB4_xy))
axis square
colorbar
title('\angle U_{MLB4}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,11)
imagesc(x,y,abs(U_Mie))
axis square
colorbar
title('|U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,12)
imagesc(x,y,angle(U_Mie))
axis square
colorbar
title('\angle U_{Mie}')
xlim(figure_window)
ylim(figure_window)

set(gcf, 'Position', get(0, 'Screensize'));


figure
subplot(3,4,1)
imagesc(x,y,abs(U_MLB_xy-U_Mie))
axis square
colorbar
title('|U_{MLB}-U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,2)
imagesc(x,y,angle(U_MLB_xy-U_Mie))
axis square
colorbar
title('\angle U_{MLB}-U_{Mie}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,3)
imagesc(x,y,abs(U_MLR_xy-U_Mie))
axis square
colorbar
title('|U_{MLR}-U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,4)
imagesc(x,y,angle(U_MLR_xy-U_Mie))
axis square
colorbar
title('\angle U_{MLR}-U_{Mie}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,7)
imagesc(x,y,abs(U_MSR_xy-U_Mie))
axis square
colorbar
title('|U_{MSR}-U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,8)
imagesc(x,y,angle(U_MSR_xy-U_Mie))
axis square
colorbar
title('\angle U_{MSR}-U_{Mie}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,5)
imagesc(x,y,abs(U_MLB2_xy-U_Mie))
axis square
colorbar
title('|U_{MLB2}-U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,6)
imagesc(x,y,angle(U_MLB2_xy-U_Mie))
axis square
colorbar
title('\angle U_{MLB2}-U_{Mie}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,9)
imagesc(x,y,abs(U_MLB4_xy-U_Mie))
axis square
colorbar
title('|U_{MLB4}-U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,10)
imagesc(x,y,angle(U_MLB4_xy-U_Mie))
axis square
colorbar
title('\angle U_{MLB4}-U_{Mie}')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,11)
imagesc(x,y,abs(U_Mie))
axis square
colorbar
title('|U_{Mie}|')
xlim(figure_window)
ylim(figure_window)

subplot(3,4,12)
imagesc(x,y,angle(U_Mie))
axis square
colorbar
title('\angle U_{Mie}')
xlim(figure_window)
ylim(figure_window)

set(gcf, 'Position', get(0, 'Screensize'));

norm(U_Mie-U_MOM,2)/norm(U_MOM,2)

figure
subplot(1,2,1)
imagesc(abs(U_Mie./U_MOM))
subplot(1,2,2)
imagesc(angle(U_Mie./U_MOM))
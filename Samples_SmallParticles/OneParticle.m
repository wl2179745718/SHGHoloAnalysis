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
thetain = 0/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
shape = 'sphere';
nsphere=1.002;%n_array(ii);%1.2;
rad = 0.1*um;
MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

V_total = VoxelCounter(shape);

%shape = 'softball';
%n_center=1.002;%1.2;
%sigma = 2*um;
%MakeSoftballHDF5(sigma, [n_center,n_imm], Box, X, Y, z,N,dr, k0);

% small box for MOM
box = [2.2*rad, 2.2*rad, 2.2*rad];

% Correct ouput wave
angle_correction = exp(1i*2*pi*(0-z(1)+dr(3))*cos(thetain)/lambda);

uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

[X1,Z1] = meshgrid(x,z);
uin_xz = exp(1i*2*pi*fxin*X1).*exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

[Y1,Z1] = meshgrid(y,z);
uin_yz = exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

x0_index = find(x==0,1);
y0_index = find(y==0,1);

dGk = -0.001;
Eps = 0;
% initialize
%udL = uin;
Pdz  = Propagator(lambda,FX,FY,dr(3));
Gamma = Gammadz(lambda,FX,FY);
Gdz  = G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);%G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
Diagnose_Green(Gdz, fx, fy, lambda);

[U_MLB_xy, U_MLB_xz, U_MLB_yz]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);

figure_window = [-Box(1)/3/um Box(1)/3/um];

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

%lmax=20;
%angle_correction = exp(1i*2*pi*(0-z(1)+dr(3))*cos(thetain)/lambda);
%U_Mie = angle_correction.*ScalarMieField(X,Y,z,thetain,lmax,n_imm,nsphere,k0,rad);
G_particle = ScalarG(X, Y, z(end), k0);

U_RemiJohn = angle_correction*V_total*G_particle;

figure
subplot(1,2,1)
imagesc(x/um,y/um,abs(U_RemiJohn))
axis square
colorbar
title('|Us_{RemiJohn}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(1,2,2)
imagesc(x/um,y/um,angle(U_RemiJohn))
axis square
colorbar
title('\angle Us_{RemiJohn}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

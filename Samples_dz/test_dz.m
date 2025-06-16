clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 1*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [40*um, 40*um, 20*um];
dr  = [ 20*nm,  20*nm,20*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
thetain = 15/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
shape = 'sphere';
nsphere=1.002;%1.2;
rad = 3*um;
MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

% Correct ouput wave
uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

[X1,Z1] = meshgrid(x,z);
uin_xz = exp(1i*2*pi*fxin*X1).*exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

[Y1,Z1] = meshgrid(y,z);
uin_yz = exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

x0_index = find(x==0,1);
y0_index = find(y==0,1);

Pdz  = Propagator(lambda,FX,FY,dr(3));
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);
Diagnose_Green(Gdz, fx, fy, lambda);
[Us_xy, Us_xz, Us_yz]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);


dielectricSphere=X.^2+Y.^2<=rad^2;
Vn=-double(dielectricSphere);
Us=myifft((myfft(uin.*Vn)).*(Gdz))*dr(1)*dr(2)*dr(3);
norm(Us,2)

dr = [dr(1),dr(2),dr(3)/2];
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
Us=myifft((myfft(uin.*Vn)).*(Gdz))*dr(1)*dr(2)*dr(3);
norm(Us,2)

% test integral
%-sum(sum(V_tst))/k0^2/(nsphere^2-1)

%4*pi/3*rad^3


% test Green's function

dielectricSphere=X.^2+Y.^2<=rad^2;
Vn=-double(dielectricSphere);
GdzVn1=myifft((myfft(uin.*Vn)).*(Gdz));

dr = [dr(1),dr(2),dr(3)/2];
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
GdzVn2=myifft((myfft(uin.*Vn)).*(Gdz));

Pdz  = Propagator(lambda,FX,FY,dr(3));
GdzVn3=myifft(Pdz.*(myfft(GdzVn2)));

errorG = norm(GdzVn1 - GdzVn3,2)/norm(GdzVn1,2)

figure_window = [-3*rad/um 3*rad/um];

figure
subplot(2,5,1)
imagesc(x/um,y/um,abs(Us_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,6)
imagesc(x/um,y/um,angle(Us_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,5,[2,7])
imagesc(z/um,x/um,abs(Us_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[3,8])
imagesc(z/um,x/um,angle(Us_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z(um)')
ylabel('x(um)')

subplot(2,5,[4,9])
imagesc(z/um,y/um,abs(Us_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z(um)')
ylabel('y(um)')

subplot(2,5,[5,10])
imagesc(z/um,y/um,angle(Us_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z(um)')
ylabel('y(um)')

set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',36)
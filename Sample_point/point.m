clear variables; close all; clc; addpath(genpath('../Functions'));Units;FigureDefault;

% Box
lambda0 = 1; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [160, 160, 5];
dr  = [0.02, 0.02, 0.01];
dz1 = 0.01;
figure_window = 20;

[Box,N] = Box_Regularization(Box,dr);

[x, y, z, x0_index, y0_index, z0_index, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

Re_fz = real(sqrt(1-FX.^2-FY.^2));
Prop_part = ones(size(Re_fz));
Prop_part(Re_fz==0)=0;

Pdz  = Propagator(lambda,FX,FY,dr(3));

z_positive = z(z0_index+1:end);
%Diagnose_Propagator_Unitary(lambda,FX,FY,dr,x0_index,fy);
%Diagnose_Propagator(X, Y, z, dr, x0_index,y0_index, z0_index, fx, fy, Pdz, k0, lambda);
%Diagnose_Propagator_accumulation(X, Y, z_positive, dr, Pdz, k0, lambda);
%Diagnose_Propagator_accumulation_clipped(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda);
%Diagnose_with_change_dz_for_the_first_layer(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda, dz1);
%Diagnose_with_change_dz1_for_the_first_layer_no_clip(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda, dz1);
%Diagnose_Greens_funtion_integral(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1);
%Diagnose_Greens_funtion_integral_NA(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1);

dz1_array = [0.01,0.1,1,10];
NA = 0.3;

%Boxsize = 40:10:320;
%Diagnose_Greens_funtion_integral_NA_2D_scan(k0, lambda, NA, Boxsize, dz1_array)

meshsize = 0.01:0.01:0.08;
Diagnose_Greens_funtion_integral_NA_2D_scan_mesh(k0, lambda, NA, meshsize, dz1_array)

%Diagnose_Propagator_accumulation_Gk(z_positive, FX, FY, fx, fy, dr, Pdz, lambda);

%Diagnose_Propagator_Accumulation_Greens_Function_Comparison(X, Y, z_positive, FX, FY, fx, fy, dr, k0, lambda);

% Generate the input plane wave
thetain = 0/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);


%Diagnose_Propagator_Plane_Wave(uin,z,dr,thetain,lambda,Pdz);

% The phaser factor of the field at the z=0 plane
angle_correction = exp(1i*2*pi*(0-z(1))*cos(thetain)/lambda);

% The sphere
shape = 'point';
nsphere=1.002;%n_array(ii);%1.2;
MakePoint([nsphere,n_imm], Box, x0_index, y0_index, z0_index, N, dr,k0);

V_total = VoxelCounter(shape);
N_Voxel = V_total./(-(k0)^2*((nsphere).^2-n_imm^2))./dr(1)./dr(2)./dr(3);

%shape = 'softball';
%n_center=1.002;%1.2;
%sigma = 2*um;
%MakeSoftballHDF5(sigma, [n_center,n_imm], Box, X, Y, z,N,dr, k0);

% small box for MOM
box = [2, 2, 2];

% Correct ouput wave

uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(end)-z(1)+dr(3))*cos(thetain)/lambda);

[X1,Z1] = meshgrid(x,z);
uin_xz = exp(1i*2*pi*fxin*X1).*exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);

[Y1,Z1] = meshgrid(y,z);
uin_yz = exp(1i*2*pi*(Z1-z(1)+dr(3))*cos(thetain)/lambda);


dGk = -0.001;
Eps = 0;
% initialize
%udL = uin;
Gamma = Gammadz(lambda,FX,FY);
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
Diagnose_Green(Gdz, fx, fy, lambda);

[U_MLB_xy, U_MLB_xz, U_MLB_yz]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);

U_MLB_xz_positive = U_MLB_xz(z0_index:end,:);

U_G_xz = 0*U_MLB_xz_positive;
parfor ii=1:size(z_positive,2)
    U_G_xz(ii,:) = angle_correction*V_total*ScalarG(X(:,x0_index), Y(:,x0_index), z_positive(ii), k0);
end


parfor ii=1:size(z_positive,2)
    error_U_xy(ii) = norm(U_G_xz(ii,:) - U_MLB_xz_positive(ii,:),2)/norm(U_G_xz(ii,:),2);
end
figure
plot(z_positive,error_U_xy)
set(gcf, 'Position', get(0, 'Screensize'));


U_P_xz = 0*U_MLB_xz_positive;
SV_src = zeros(N(1),N(2));
SV_src(x0_index,y0_index) = angle_correction*V_total;

Us=myifft((myfft(SV_src)).*Pdz);

for ii = 1:size(z_positive,2)
    U_P_xz(ii,:) = Us(:,x0_index);
    Us=myifft((myfft(Us)).*Pdz);
    error_UP_xy(ii) = norm(U_P_xz(ii,:) - U_MLB_xz_positive(ii,:),2)/norm(U_P_xz(ii,:),2);
end
figure
plot(z_positive,error_UP_xy)
set(gcf, 'Position', get(0, 'Screensize'));






U_P_xz = 0*U_MLB_xz_positive;

SV_src = zeros(N(1),N(2));
SV_src(x0_index,y0_index) = angle_correction*V_total;
for ii = 1:size(z_positive,2)
    Prop = Propagator(lambda,FX,FY,z(ii));
    Us=myifft((myfft(SV_src)).*Prop);
    U_P_xz(ii,:) = Us(:,x0_index);
    error_UP_xy(ii) = norm(U_P_xz(ii,:) - U_MLB_xz_positive(ii,:),2)/norm(U_P_xz(ii,:),2);
end

figure('Name','xz plane')
plot(z_positive,error_UP_xy)
set(gcf, 'Position', get(0, 'Screensize'));
z_positive = z(z0_index+1:end);
figure
subplot(1,3,1)
imagesc(z_positive,y,abs(U_MLB_xz_positive)')
%axis square
colorbar
title('MLB')

subplot(1,3,2)
imagesc(z_positive,y,abs(U_G_xz)')
%axis square
colorbar
title('Greens function')

subplot(1,3,3)
imagesc(z_positive,y,abs(U_MLB_xz_positive-U_G_xz)')
%axis square
colorbar
title('Difference')
set(gcf, 'Position', get(0, 'Screensize'));


figure('Name','log field')
subplot(1,3,1)
imagesc(z_positive,y,log( abs(U_MLB_xz_positive)' ))
ylim([-figure_window/2 figure_window/2])
axis equal
colorbar
title('MLB')

subplot(1,3,2)
imagesc(z_positive,y,log( abs(U_G_xz)' ))
ylim([-figure_window/2 figure_window/2])
axis equal
colorbar
title('Greens function')

subplot(1,3,3)
imagesc(z_positive,y,log( abs(U_MLB_xz_positive-U_G_xz)' ))
ylim([-figure_window/2 figure_window/2])
axis equal
colorbar
title('Difference')
set(gcf, 'Position', get(0, 'Screensize'));



for ii = 1:size(z_positive,2)
    ii
    errorU(ii) = norm( U_G_xz(z0_index+1:end,:)-U_P_xz(z0_index+1:end,:),2 ) / norm( U_G_xz(z0_index+1:end,:),2 );
end

figure
plot(z_positive,errorU)
set(gcf, 'Position', get(0, 'Screensize'));

figure_window = [-Box(1)/3 Box(1)/3];

figure
subplot(2,5,1)
imagesc(x,y,abs(U_MLB_xy))
axis square
colorbar
title('|Us_{out}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

subplot(2,5,6)
imagesc(x,y,angle(U_MLB_xy))
axis square
colorbar
title('\angle Us_{out}')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

subplot(2,5,[2,7])
imagesc(z,x,abs(U_MLB_xz.'))
axis square
colorbar
title('|Us_{xz}|')
axis equal
xlabel('z')
ylabel('x')

subplot(2,5,[3,8])
imagesc(z,x,angle(U_MLB_xz.'))
axis square
colorbar
title('\angle Us_{xz}')
axis equal
xlabel('z')
ylabel('x')

subplot(2,5,[4,9])
imagesc(z,y,abs(U_MLB_yz.'))
axis square
colorbar
title('|Us_{yz}|')
axis equal
xlabel('z')
ylabel('y')

subplot(2,5,[5,10])
imagesc(z,y,angle(U_MLB_yz.'))
axis square
colorbar
title('\angle Us_{yz}')
axis equal
xlabel('z')
ylabel('y')

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
imagesc(x,y,abs(U_RemiJohn))
axis square
colorbar
title('|Us_{RemiJohn}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

subplot(1,2,2)
imagesc(x,y,angle(U_RemiJohn))
axis square
colorbar
title('\angle Us_{RemiJohn}')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

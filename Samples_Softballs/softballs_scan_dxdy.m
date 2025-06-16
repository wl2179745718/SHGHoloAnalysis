clear variables; close all; clc; addpath(genpath('../Functions'));Units;FigureDefault;

% Box
lambda0 = 1; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

NA = 0.99;

dxy_array = 0.02:0.004:0.06;

for jj = 1:size(dxy_array,2)
    dxy_array(jj)

Box = [80, 80, 10];
dr  = [dxy_array(jj), dxy_array(jj), 0.04];
dz1 = 0.01;
figure_window = 20;

[Box,N] = Box_Regularization(Box,dr);

[x, y, z, x0_index, y0_index, z0_index, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

Re_fz = real(sqrt(1-FX.^2-FY.^2));
Prop_part = ones(size(Re_fz));
Prop_part(Re_fz==0)=0;

Pdz  = Propagator(lambda,FX,FY,dr(3));

%z_positive = z(z0_index+1:end);
%Diagnose_Propagator_Unitary(lambda,FX,FY,dr,x0_index,fy);
%Diagnose_Propagator(X, Y, z, dr, x0_index,y0_index, z0_index, fx, fy, Pdz, k0, lambda);
%Diagnose_Propagator_accumulation(X, Y, z_positive, dr, Pdz, k0, lambda);
%Diagnose_Propagator_accumulation_clipped(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda);
%Diagnose_with_change_dz_for_the_first_layer(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda, dz1);
%Diagnose_with_change_dz1_for_the_first_layer_no_clip(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda, dz1);
%Diagnose_Greens_funtion_integral(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1);
%Diagnose_Greens_funtion_integral_NA(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1);

%dz1_array = [0.01,0.1,1,10];
%NA = 0.3;

%Boxsize = 40:10:320;
%Diagnose_Greens_funtion_integral_NA_2D_scan(k0, lambda, NA, Boxsize, dz1_array)

%meshsize = 0.01:0.01:0.08;
%Diagnose_Greens_funtion_integral_NA_2D_scan_mesh(k0, lambda, NA, meshsize, dz1_array)

%Diagnose_Propagator_accumulation_Gk(z_positive, FX, FY, fx, fy, dr, Pdz, lambda);

%Diagnose_Propagator_Accumulation_Greens_Function_Comparison(X, Y, z_positive, FX, FY, fx, fy, dr, k0, lambda);

% Generate the input plane wave
thetain = 0/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);


%Diagnose_Propagator_Plane_Wave(uin,z,dr,thetain,lambda,Pdz);

% The phaser factor of the field at the z=0 plane
%angle_correction = exp(1i*2*pi*(0-z(1))*cos(thetain)/lambda);

% The sphere
%shape = 'point';
%nsphere=1.002;%n_array(ii);%1.2;
%MakePoint([nsphere,n_imm], Box, x0_index, y0_index, z0_index, N, dr,k0);

%shape = 'softball';
%n_center=1.002;%1.2;
%sigma = 2*um;
%MakeSoftballHDF5(sigma, [n_center,n_imm], Box, X, Y, z,N,dr, k0);

shape = 'softballs';
rj = [0 0 0; 0 1 0; 1 0 0].';
sigma0 = 1/3;
sigma = sigma0*[1 1 1];
n0 = 1.01;
n = n0*[1 1 1];
MakeSoftballArray(rj, sigma, n,n_imm, Box, X, Y, z,N,dr, k0);



%V_total = VoxelCounter(shape);
%N_Voxel = V_total./(-(k0)^2*((nsphere).^2-n_imm^2))./dr(1)./dr(2)./dr(3);

% small box for MOM
%box = [2, 2, 2];

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



Gdz  = G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);%G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
%Diagnose_Green(Gdz, fx, fy, lambda);

[U_MLB_xy, U_MLB_xz, U_MLB_yz]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);

U_MLB_xz_positive = U_MLB_xz(z0_index:end,:);


KX = 2*pi*FX;
KY = 2*pi*FY;

PZ = Propagator(lambda,FX,FY, -z(end));
E_MLB_image=myifft(PZ.*(myfft(U_MLB_xy)));

%E_MLB_image = BPM(U_MLB_xy, z(end), NA, k0, KX, KY);

Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(FX,FY,n_imm,lambda,  dr(3),dGk,Eps)./dr(1)./dr(2);GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);
%Diagnose_Green(Gdz, fx, fy, lambda);

[U_MLB_xy_Oli, U_MLB_xz_Oli, U_MLB_yz_Oli]=MLB(shape,uin,Gdz,Pdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index);

U_MLB_xz_positive_Oli = U_MLB_xz_Oli(z0_index:end,:);


E_MLB_image_Oli=myifft(PZ.*(myfft(U_MLB_xy_Oli)));

%E_MLB_image = BPM(U_MLB_xy, z(end), NA, k0, KX, KY);

ki = 2*pi*[0 0 1];
KZ = sqrt(k0^2-KX.^2-KY.^2);
Gaus_F = (n0^2 - 1)*(2*pi)^1.5*exp(- sigma(1)^2/2 *  ( (KX-ki(1)).^2 + (KY-ki(2)).^2 + (KZ-ki(3)).^2 ) );
Gaus_F(FX.^2 + FY.^2 > 0.99^2) = 0;
Ek_single = 1i*k0^2/4/pi./sqrt(k0^2-KX.^2-KY.^2+1e-12i).* Gaus_F;
%Ek_single = 1i*k0^2/4/pi./sqrt(k0^2-KX.^2-KY.^2+1e-2i).* Gaus_F;

ArrayF = 0;
for ii = 1:size(rj, 2)
ArrayF = ArrayF + exp(1i*(rj(1,ii)*(ki(1)-KX)+rj(2,ii)*(ki(2)-KY)+rj(3,ii)*(ki(3)-KZ)  ));

%1+exp(1i*(rj(1,2)*(ki(1)-KX)+rj(2,2)*(ki(2)-KY)+rj(3,2)*(ki(3)-KZ)  ));
end

Ek_all = Ek_single.*ArrayF;

Er = 1/dr(1)/dr(2)*sigma0^3*2*pi*myifft(Ek_all);


figure_window = [-5 5];



figure('Name','xz and yz field')
subplot(2,3,1)
imagesc(x,y,abs(E_MLB_image_Oli))
axis square
colorbar
title('|E_{Oli}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

subplot(2,3,4)
imagesc(x,y,wrapToPi(angle(E_MLB_image_Oli)))
axis square
colorbar
title('\angle E_{Oli}')
xlim(figure_window)
ylim(figure_window)
clim([-pi pi])
xlabel('x')
ylabel('y')

subplot(2,3,2)
imagesc(x,y,abs(E_MLB_image))
axis square
colorbar
title('|E_{MLB}|')
xlim(figure_window)
ylim(figure_window)
xlabel('x')
ylabel('y')

subplot(2,3,5)
imagesc(x,y,wrapToPi(angle(E_MLB_image)))
axis square
colorbar
title('\angle E_{MLB}')
xlim(figure_window)
ylim(figure_window)
clim([-pi pi])
xlabel('x')
ylabel('y')

subplot(2,3,3)
imagesc(x,y,abs(Er))
xlim(figure_window)
ylim(figure_window)
axis square
colorbar
title('|E_{k}|')
xlabel('x')
ylabel('y')

subplot(2,3,6)
imagesc(x,y,wrapToPi(angle(Er)))
xlim(figure_window)
ylim(figure_window)
clim([-pi pi])
axis square
colorbar
title('\angle E_{k}')
xlabel('x')
ylabel('y')

sgtitle('E(r_\perp,z=0)')
set(gcf, 'Position', get(0, 'Screensize'));


% power error
power_error_scan_NA(fx, FX, FY, dfx, dfy, Re_fz, E_MLB_image_Oli, E_MLB_image, Er)

filter_out = abs( sqrt(FX.^2+FY.^2) ) < NA;
Er_power = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(myfft(Er)).^2 ));
error_power_Oli = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(myfft(Er - E_MLB_image_Oli)).^2 ));
error_power_MLB = dfx * dfy * sum(sum( filter_out.*Re_fz.*abs(myfft(Er - E_MLB_image)).^2 ));
err_pwr_ratio_Oli(jj) = error_power_Oli/Er_power;
err_pwr_ratio_MLB(jj) = error_power_MLB/Er_power;

end

figure
hold on
plot(dxy_array, err_pwr_ratio_Oli)
plot(dxy_array, err_pwr_ratio_MLB)
xlabel('dx = dy')
ylabel('error')
legend('Olivier','MLB')

sgtitle('E(r_\perp,z=0)')
set(gcf, 'Position', get(0, 'Screensize'));
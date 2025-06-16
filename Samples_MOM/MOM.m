clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 10*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [300*um, 300*um, 50*um];
dr  = [ 708*nm,  708*nm,128*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

% Generate the input plane wave
thetain = 15/180*pi; % input angle
[thetain,fxin] = Angle_Regularization(thetain,lambda,dfx);
uin = exp(1i*2*pi*fxin*X);

% The sphere
shape = 'sphere';
nsphere=1.02;%1.2;
rad = 10*um;
MakeSphereHDF5(rad, [nsphere,n_imm], Box, X, Y, z, N, dr,k0);

% small box
box = [2.1*rad, 2.1*rad, 2.1*rad];

box_xy_index = find(abs(X)<box(1)/2 & abs(Y)<box(2)/2);
box_X = X(box_xy_index);
box_Y = Y(box_xy_index);
box_z_index  = find(abs(z)<box(3)/2);
box_z = z(box_z_index);

% MOM
chunk_size = [N(1),  N(2)];

U_MOM_1d = zeros(1,N(1)*N(2));

for ii = 1:size(box_z_index,2)
    ii
    % Read Vn
    start=[1 1 box_z_index(ii)]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read('../Medium/sphere.h5','/../Medium/sphere',start,count);
    Vn = Vn(box_xy_index);
    % Calculate incident field
    Uin_n = exp(1i*2*pi*fxin*box_X)*exp(1i*2*pi*(box_z(ii)-z(1)+dr(3))*cos(thetain)/lambda);
    % Calculate the Green's function
    parfor jj = 1:N(1)*N(2)
        Rn = sqrt( (X(jj)-box_X).^2 + (Y(jj)-box_Y).^2 + (z(end)-box_z(ii))^2  );
        Gn = -1/4/pi*exp(1i*k0*Rn)./Rn;
        U_MOM_1d(jj) = U_MOM_1d(jj) + sum( Vn.*Uin_n.*Gn );
    end
end

U_MOM = dr(1)*dr(2)*dr(3)*reshape(U_MOM_1d, N(1), N(2));

figure
subplot(1,2,1)
imagesc(x,y,abs(U_MOM))
axis square
colorbar
title('|U_{MOM}|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(1,2,2)
imagesc(x,y,angle(U_MOM))
axis square
colorbar
title('\angle U_{MOM}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
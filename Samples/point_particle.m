clear variables; close all; clc; addpath(genpath('../Functions'));

ps          = 0.08;                 % pixel size (x,y,z) in object space (microns)
lambda      = 1;                  % central wavelength (microns)  
n_imm       = 1;                % refractive index of immersion media
nsphere     = 1.02;%1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
L = [20 20 20];
d = [1 0 0];
%N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
%L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
rad = 1.5*ps;
defocus = [0];

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.25;

U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

ord = 1;

E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MSL_3d=MultiSlice(U_inp,RI,lambda,ps,fxx,fyy,n_imm,'Vol');

V_inc = 0*V;
RI_inc = n_imm * ones(size(RI));
E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
E_MSL_3d_inc=MultiSlice(U_inp,RI_inc,lambda,ps,fxx,fyy,n_imm,'Vol');

E_MSR_3d_sca= E_MSR_3d - E_MSR_3d_inc;
E_MLR_3d_sca= E_MLR_3d - E_MLR_3d_inc;
E_MLB_3d_sca= E_MLB_3d - E_MLB_3d_inc;
E_MSL_3d_sca= E_MSL_3d - E_MSL_3d_inc;

E_MSR_3d_sca_xy = squeeze(E_MSR_3d_sca(:,:,end));
E_MLR_3d_sca_xy = squeeze(E_MLR_3d_sca(:,:,end));
E_MLB_3d_sca_xy = squeeze(E_MLB_3d_sca(:,:,end));
E_MSL_3d_sca_xy = squeeze(E_MSL_3d_sca(:,:,end));

E_MSR_3d_sca_yz = squeeze(E_MSR_3d_sca(end,:,:));
E_MLR_3d_sca_yz = squeeze(E_MLR_3d_sca(end,:,:));
E_MLB_3d_sca_yz = squeeze(E_MLB_3d_sca(end,:,:));
E_MSL_3d_sca_yz = squeeze(E_MSL_3d_sca(end,:,:));

E_MSR_3d_sca_xz = squeeze(E_MSR_3d_sca(:,end,:));
E_MLR_3d_sca_xz = squeeze(E_MLR_3d_sca(:,end,:));
E_MLB_3d_sca_xz = squeeze(E_MLB_3d_sca(:,end,:));
E_MSL_3d_sca_xz = squeeze(E_MSL_3d_sca(:,end,:));

% xy plane

[theta1D_xy, R_xy] = XYz2thetaR(d, X, Y, z(end));

c=299792458;
[~,~,~,ETheta1D_xy] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,theta1D_xy);
ETheta_xy = reshape(ETheta1D_xy, [size(X,1) size(X,2)]);
E_vec_Mie = ETheta_xy./( exp(1i*k*z(end))/z(end) ).*( exp(1i*k*R_xy)./R_xy );


l_max = 40;
[E_sca_Mie_1D] = scalarMie(l_max,k,nsphere,rad,theta1D_xy);
E_sca_Mie = reshape(E_sca_Mie_1D, size(X));

figure
subplot(2,2,1)
imagesc(x,y,abs(E_sca_Mie))
axis square
colorbar
title('|E_{scaMie}|')

subplot(2,2,2)
imagesc(x,y,abs(E_vec_Mie))
axis square
colorbar
title('|E_{vecMie}|')

subplot(2,2,3)
imagesc(x,y,angle(E_sca_Mie))
axis square
colorbar
title('\angle E_{scaMie}')

subplot(2,2,4)
imagesc(x,y,angle(E_vec_Mie))
axis square
colorbar
title('\angle E_{vecMie}')

f_title = sprintf('Mie xy plane');
sgtitle(f_title)

figure
subplot(2,4,1)
imagesc(x,y,abs(E_MSR_3d_sca_xy))
axis square
colorbar
title('|E_{MSR}|')
subplot(2,4,5)
imagesc(x,y,angle(E_MSR_3d_sca_xy))
axis square
colorbar
title('\angle E_{MSR}')

subplot(2,4,2)
imagesc(x,y,abs(E_MLR_3d_sca_xy))
axis square
colorbar
title('|E_{MLR}|')
subplot(2,4,6)
imagesc(x,y,angle(E_MLR_3d_sca_xy))
axis square
colorbar
title('\angle E_{MLR}')

subplot(2,4,3)
imagesc(x,y,abs(E_MLB_3d_sca_xy))
axis square
colorbar
title('|E_{MLB}|')
subplot(2,4,7)
imagesc(x,y,angle(E_MLB_3d_sca_xy))
axis square
colorbar
title('\angle E_{MLB}')

subplot(2,4,4)
imagesc(x,y,abs(E_MSL_3d_sca_xy))
axis square
colorbar
title('|E_{MSL}|')
subplot(2,4,8)
imagesc(x,y,angle(E_MSL_3d_sca_xy))
axis square
colorbar
title('\angle E_{MSL}')

f_title = sprintf('Simu xy plane');
sgtitle(f_title)

theta_xy = reshape(theta1D_xy,size(X));
figure
imagesc(x,y,theta_xy)
axis square
colorbar


X_1D = reshape(X,size(theta1D_xy));
Y_1D = reshape(Y,size(theta1D_xy));
for ii = 1:size(X_1D,2)
    E = dyadic_Green([X_1D(ii);Y_1D(ii);z(end)],k,d);
    E_dyGreen_1D(ii) = norm(E);
end

E_dyGreen = reshape(E_dyGreen_1D,size(X));

E_scGreen = scalar_Green(k,R_xy);

figure
subplot(1,3,1)
imagesc(x,y,E_dyGreen)
axis square
colorbar
title('|E_{dyGreen}|_{xy}')
subplot(1,3,2)
imagesc(x,y,abs(E_scGreen))
axis square
colorbar
title('|E_{scGreen}|_{xy}')
subplot(1,3,3)
imagesc(x,y,angle(E_scGreen))
axis square
colorbar
title('\angle E_{scGreen}_{xy}')

f_title = sprintf('Green xy');
sgtitle(f_title)

% xz plane

[X,Z]=meshgrid(x,z);

[theta1D_xz, R_xz] = XZy2thetaR(d, X, y(end), Z);

c=299792458;
[~,~,~,ETheta1D_xz] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,theta1D_xz);
ETheta_xz = reshape(ETheta1D_xz, [size(X,1) size(X,2)]);
E_vec_Mie = ETheta_xz./( exp(1i*k*z(end))/z(end) ).*( exp(1i*k*R_xz)./R_xz );


l_max = 40;
[E_sca_Mie_1D] = scalarMie(l_max,k,nsphere,rad,theta1D_xz);
E_sca_Mie = reshape(E_sca_Mie_1D, size(X));

figure
subplot(2,2,1)
imagesc(x,z,abs(E_sca_Mie))
axis square
colorbar
title('|E_{scaMie}|')

subplot(2,2,2)
imagesc(x,z,abs(E_vec_Mie))
axis square
colorbar
title('|E_{vecMie}|')

subplot(2,2,3)
imagesc(x,z,angle(E_sca_Mie))
axis square
colorbar
title('\angle E_{scaMie}')

subplot(2,2,4)
imagesc(x,z,angle(E_vec_Mie))
axis square
colorbar
title('\angle E_{vecMie}')

f_title = sprintf('Mie xz plane');
sgtitle(f_title)

figure
subplot(2,4,1)
imagesc(x,z,abs(E_MSR_3d_sca_xz))
axis square
colorbar
title('|E_{MSR}|')
subplot(2,4,5)
imagesc(x,z,angle(E_MSR_3d_sca_xz))
axis square
colorbar
title('\angle E_{MSR}')

subplot(2,4,2)
imagesc(x,z,abs(E_MLR_3d_sca_xz))
axis square
colorbar
title('|E_{MLR}|')
subplot(2,4,6)
imagesc(x,z,angle(E_MLR_3d_sca_xz))
axis square
colorbar
title('\angle E_{MLR}')

subplot(2,4,3)
imagesc(x,z,abs(E_MLB_3d_sca_xz))
axis square
colorbar
title('|E_{MLB}|')
subplot(2,4,7)
imagesc(x,z,angle(E_MLB_3d_sca_xz))
axis square
colorbar
title('\angle E_{MLB}')

subplot(2,4,4)
imagesc(x,z,abs(E_MSL_3d_sca_xz))
axis square
colorbar
title('|E_{MSL}|')
subplot(2,4,8)
imagesc(x,z,angle(E_MSL_3d_sca_xz))
axis square
colorbar
title('\angle E_{MSL}')

f_title = sprintf('Simu xz plane');
sgtitle(f_title)

theta_xz = reshape(theta1D_xz,size(X));
figure
imagesc(x,z,theta_xz)
axis square
colorbar

X_1D = reshape(X,size(theta1D_xy));
Z_1D = reshape(Z,size(theta1D_xy));
for ii = 1:size(X_1D,2)
    E = dyadic_Green([X_1D(ii);y(end);Z_1D(ii)],k,d);
    E_dyGreen_1D(ii) = norm(E);
end

E_dyGreen = reshape(E_dyGreen_1D,size(X));

E_scGreen = scalar_Green(k,R_xz);

figure
subplot(1,3,1)
imagesc(x,z,E_dyGreen)
axis square
colorbar
title('|E_{dyGreen}|_{xz}')
subplot(1,3,2)
imagesc(x,z,abs(E_scGreen))
axis square
colorbar
title('|E_{scGreen}|_{xz}')
subplot(1,3,3)
imagesc(x,z,angle(E_scGreen))
axis square
colorbar
title('\angle E_{scGreen}_{xz}')

f_title = sprintf('Green xz');
sgtitle(f_title)









clear variables; close all; clc; addpath(genpath('../Functions'));

% Global Parameters

dx          = .1;       % pixel size (x,y) in object space (microns)
dz          = .1;       % pixel size (x,y) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
%N           = [N_array(ii)-1, N_array(ii)-1, 2^9];                  % lateral pixel dimension 
L = [30, 30, 30];
delta = [dx, dx, dz];
%N = round(L./delta);
deltaf = 1./L;     

% The grid

[N, x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

% The sphere

shape = 'sphere';

nsphere=1.01;%1.2;
n=[nsphere,n_imm];
rad = 2;
MakeSphereHDF5(rad, n, L, X, Y, z, N, delta,k0);

% The incident field

NA_in = 0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);  
U_in=exp(1i*k*sin(phi)*X);

% The Green's function

dGk = 1;
Eps=n_imm^2/lambda^2*0.05;
Gdz = G_kx_ky(fxx,fyy,n_imm,lambda,dz,dGk,Eps);
Pdz = Propagator(n_imm,lambda,fxx,fyy,dz);
Gamma = Gammadz(n_imm,lambda,fxx,fyy);

order = 1;

E_MLB_3d=MultiLayerBornv3(shape,dz,0,U_in,Gdz,Pdz);
E_MLR_3d=MultiLayerRytovv3(shape,dz,0,U_in,Gdz,Pdz);
E_MSR_3d=MultiSlabRytovv3(shape,dz,0,U_in,order,Gdz,Pdz,Gamma);
E_MSL_3d=MultiSliceV3(shape,U_in,Pdz,0,k0,dz,n_imm);

E_MLB_3d_inc=MultiLayerBornv3(shape,dz,1,U_in,Gdz,Pdz);
E_MLR_3d_inc=MultiLayerRytovv3(shape,dz,1,U_in,Gdz,Pdz);
E_MSR_3d_inc=MultiSlabRytovv3(shape,dz,1,U_in,order,Gdz,Pdz,Gamma);
E_MSL_3d_inc=MultiSliceV3(shape,U_in,Pdz,1,k0,dz,n_imm);

E_sca_MSR = E_MSR_3d-E_MSR_3d_inc;
E_sca_MLR = E_MLR_3d-E_MLR_3d_inc;
E_sca_MLB = E_MLB_3d-E_MLB_3d_inc;
E_sca_MSL = E_MSL_3d-E_MSL_3d_inc;

E_sca_MLB = E_sca_MLB./max(max(E_sca_MLB));
E_sca_MLR = E_sca_MLR./max(max(E_sca_MLR));
E_sca_MSR = E_sca_MSR./max(max(E_sca_MSR));
E_sca_MSL = E_sca_MSL./max(max(E_sca_MSL));


[E_sca_Mie] = Mie_plane(X,Y,z(end),phi,k,rad,c/lambda,n_imm^2,1,nsphere^2,1,40);


figure
subplot(2,2,1)
hold on
plot(unwrap(angle(E_sca_MLB(round(end/2),:))))
plot(unwrap(angle(E_sca_Mie(round(end/2),:))))
subplot(2,2,2)
plot(unwrap(angle(E_sca_MLB(:,round(end/2)))))
plot(unwrap(angle(E_sca_Mie(:,round(end/2)))),'o')
subplot(2,2,3)
plot(abs(E_sca_MLB(:,round(end/2))))
plot(abs(E_sca_Mie(:,round(end/2))),'o')
subplot(2,2,4)
plot(abs(E_sca_MLB(round(end/2),:)))
plot(abs(E_sca_Mie(round(end/2),:)),'o')

%[~,m_index]=max(abs(E_sca_MLB(round(end/2),:)));
%E_sca_MLB(round(end/2),m_index)
%x(m_index)

%[m_value,m_index]=max(abs(E_sca_Mie(round(end/2),:)));
%E_sca_Mie(round(end/2),m_index)
%x(m_index)

%corr_m = corr(E_sca_MLB(round(end/2),:), E_sca_Mie(round(end/2),:));
%crr = xcorr2(E_sca_MLB,E_sca_Mie);
%[ssr,snd] = max(crr(:));
%[i,j] = ind2sub(size(crr),snd);


thrsh = 0.2;

%cond = abs(E_sca_Mie) > thrsh & abs(E_sca_MSL) > thrsh;

%E_Mie_corr = E_sca_Mie.*cond;
%E_MSL_corr = E_sca_MSL.*cond;

%fun = @(x)norm(abs((x(1)+1i*x(2))*E_Mie_corr - E_MSL_corr));
%x0 = [1,0];
%x_factor = fminsearch(fun,x0);

%E_Mie_MSL = (x_factor(1)+1i*x_factor(2))*E_sca_Mie;

E_sca_MLB = normalize(E_sca_MLB, E_sca_Mie, thrsh);
E_sca_MLR = normalize(E_sca_MLR, E_sca_Mie, thrsh);
E_sca_MSR = normalize(E_sca_MSR, E_sca_Mie, thrsh);
E_sca_MSL = normalize(E_sca_MSL, E_sca_Mie, thrsh);


NA = 0.2790;

E_MLB_sca = BPM(E_sca_MLB, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MLR_sca = BPM(E_sca_MLR, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR_sca = BPM(E_sca_MSR, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSL_sca = BPM(E_sca_MSL, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_Mie_sca = BPM(E_sca_Mie, 0, NA, k, 2*pi*fxx, 2*pi*fyy);

E_MLB_tot = BPM(E_MLB_3d, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MLR_tot = BPM(E_MLR_3d, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR_tot = BPM(E_MSR_3d, 0, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSL_tot = BPM(E_MSL_3d, 0, NA, k, 2*pi*fxx, 2*pi*fyy);


cmax_sca_E = max(max(abs(E_MLB_sca)));
cmax_sca_I = max(max(abs(E_MLB_sca).^2));
cmax_tot_E = max(max(abs(E_MLB_tot)));
cmax_tot_I = max(max(abs(E_MLB_tot).^2));
cmin_tot_E = min(min(abs(E_MLB_tot)));
cmin_tot_I = min(min(abs(E_MLB_tot).^2));

figure
%set(gcf, 'Position', get(0, 'Screensize'));
colormap(gray(256));
subplot(5, 6, 1)
imagesc(x,y,abs(E_MLB_sca))
axis square
colorbar
title('MLB |E| from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_E]);

subplot(5, 6, 2)
imagesc(x,y,angle(E_MLB_sca))
axis square
colorbar
title('MLB \angle E from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 3)
imagesc(x,y,abs(E_MLB_sca).^2)
axis square
colorbar
title('MLB I from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_I]);

subplot(5, 6, 4)
imagesc(x,y,abs(E_MLB_tot))
axis square
colorbar
title('MLB |E| from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_E cmax_tot_E]);

subplot(5, 6, 5)
imagesc(x,y,angle(E_MLB_tot))
axis square
colorbar
title('MLB \angle E from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 6)
imagesc(x,y,abs(E_MLB_tot).^2)
axis square
colorbar
title('MLB I from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_I cmax_tot_I]);

subplot(5, 6, 7)
imagesc(x,y,abs(E_MLR_sca))
axis square
colorbar
title('MLR |E| from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_E]);

subplot(5, 6, 8)
imagesc(x,y,angle(E_MLR_sca))
axis square
colorbar
title('MLR \angle E from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 9)
imagesc(x,y,abs(E_MLR_sca).^2)
axis square
colorbar
title('MLR I from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_I]);

subplot(5, 6, 10)
imagesc(x,y,abs(E_MLR_tot))
axis square
colorbar
title('MLR |E| from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_E cmax_tot_E]);

subplot(5, 6, 11)
imagesc(x,y,angle(E_MLR_tot))
axis square
colorbar
title('MLR \angle E from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 12)
imagesc(x,y,abs(E_MLR_tot).^2)
axis square
colorbar
title('MLR I from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_I cmax_tot_I]);

subplot(5, 6, 13)
imagesc(x,y,abs(E_MSR_sca))
axis square
colorbar
title('MSR |E| from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_E]);

subplot(5, 6, 14)
imagesc(x,y,angle(E_MSR_sca))
axis square
colorbar
title('MSR \angle E from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 15)
imagesc(x,y,abs(E_MSR_sca).^2)
axis square
colorbar
title('MSR I from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_I]);

subplot(5, 6, 16)
imagesc(x,y,abs(E_MSR_tot))
axis square
colorbar
title('MSR |E| from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_E cmax_tot_E]);

subplot(5, 6, 17)
imagesc(x,y,angle(E_MSR_tot))
axis square
colorbar
title('MSR \angle E from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 18)
imagesc(x,y,abs(E_MSR_tot).^2)
axis square
colorbar
title('MSR I from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_I cmax_tot_I]);

subplot(5, 6, 19)
imagesc(x,y,abs(E_MSL_sca))
axis square
colorbar
title('MSL |E| from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_E]);

subplot(5, 6, 20)
imagesc(x,y,angle(E_MSL_sca))
axis square
colorbar
title('MSL \angle E from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 21)
imagesc(x,y,abs(E_MSL_sca).^2)
axis square
colorbar
title('MSL I from E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([0 cmax_sca_I]);

subplot(5, 6, 22)
imagesc(x,y,abs(E_MSL_tot))
axis square
colorbar
title('MSL |E| from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_E cmax_tot_E]);

subplot(5, 6, 23)
imagesc(x,y,angle(E_MSL_tot))
axis square
colorbar
title('MSL \angle E from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 24)
imagesc(x,y,abs(E_MSL_tot).^2)
axis square
colorbar
title('MSL I from E_{tot}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])
clim([cmin_tot_I cmax_tot_I]);

subplot(5, 6, 25)
imagesc(x,y,abs(E_Mie_sca))
axis square
colorbar
title('Mie |E|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 26)
imagesc(x,y,angle(E_Mie_sca))
axis square
colorbar
title('Mie \angle E')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 6, 27)
imagesc(x,y,abs(E_Mie_sca).^2)
axis square
colorbar
title('Mie I')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])


angle_MLB = unwrap(angle(E_MLB_sca(round(end/2),:)));
angle_MLR = unwrap(angle(E_MLR_sca(round(end/2),:)));
angle_MSR = unwrap(angle(E_MSR_sca(round(end/2),:)));
angle_Mie = unwrap(angle(E_Mie_sca(round(end/2),:)));
% remornalize
angle_MLB = angle_MLB - angle_MLB(round(end/2));
angle_MLR = angle_MLR - angle_MLR(round(end/2));
angle_MSR = angle_MSR - angle_MSR(round(end/2));
angle_Mie = angle_Mie - angle_Mie(round(end/2));

subplot(5, 6, [28 29 30])
hold on
plot(x, angle_MLB,'LineWidth',4);
plot(x, angle_MLR,'LineStyle',"--",'LineWidth',4);
plot(x, angle_MSR,'LineStyle',":",'LineWidth',4);
plot(x, angle_Mie,'LineStyle',"-.",'LineWidth',4);
legend('MLB','MLR','MSR','Mie')
xlim([-2*rad 2*rad])

saveas(gcf,'Report20240607.png')

disp('job done')

%saveas(gcf,'Figures/NA0p0.png')
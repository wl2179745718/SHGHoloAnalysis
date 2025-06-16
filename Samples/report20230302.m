clear variables; close all; clc; addpath(genpath('../Functions'));
% unit = m
ps          = .105/2.*1e-6;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.7.*1e-6;                  % central wavelength (microns)  
n_imm       = 1.5678;                % refractive index of immersion media
nsphere=1.5717;%1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.4;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
rad = 15/2.*1e-6;
defocus = [0 1 2].*1e-6;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.25;

U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

ord = 1;

E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MSL_3d=MultiSlice(U_inp,RI,lambda,ps,fxx,fyy,n_imm);

V_inc = 0*V;
RI_inc = n_imm * ones(size(RI));
E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
E_MSL_3d_inc=MultiSlice(U_inp,RI_inc,lambda,ps,fxx,fyy,n_imm);
U_inp_end_MSR = E_MSR_3d_inc(:,:,end);
U_inp_end_MLR = E_MLR_3d_inc(:,:,end);
U_inp_end_MLB = E_MLB_3d_inc(:,:,end);
U_inp_end_MSL = E_MSL_3d_inc(:,:,end);


E_MSR_xz = squeeze(E_MSR_3d(end/2,:,:));
E_MLR_xz = squeeze(E_MLR_3d(end/2,:,:));
E_MLB_xz = squeeze(E_MLB_3d(end/2,:,:));
E_MSL_xz = squeeze(E_MSL_3d(end/2,:,:));


[E_tot_MSR,E_sca_MSR] = tot2sca(E_MSR_3d,U_inp_end_MSR);
[E_tot_MLR,E_sca_MLR] = tot2sca(E_MLR_3d,U_inp_end_MLR);
[E_tot_MLB,E_sca_MLB] = tot2sca(E_MLB_3d,U_inp_end_MLB);
[E_tot_MSL,E_sca_MSL] = tot2sca(E_MSL_3d,U_inp_end_MSL);

%E_tot_MSR = E_tot_MSR./E_tot_MSR(end/2,end/2);
%E_tot_MLR = E_tot_MLR./E_tot_MLR(end/2,end/2);
%E_tot_MLB = E_tot_MLB./E_tot_MLB(end/2,end/2);

%E_sca_MSR = E_sca_MSR./E_sca_MSR(end/2,end/2);
%E_sca_MLR = E_sca_MLR./E_sca_MLR(end/2,end/2);
%E_sca_MLB = E_sca_MLB./E_sca_MLB(end/2,end/2);

%theta = atan( sqrt(X.^2 + Y.^2)./z(end) );
%theta1D = reshape(theta, [1 size(theta,1)*size(theta,2)]);
X_1d = reshape(X, [1 size(X,1)*size(X,2)]);
Y_1d = reshape(Y, [1 size(Y,1)*size(Y,2)]);
theta1D = 0*X_1d;
V1 = [tan(phi) 0 1];
for ii = 1:size(theta1D,2)
    V2 = [X_1d(ii) Y_1d(ii) z(end)];
    theta1D(ii) = atan2(norm(cross(V1,V2)),dot(V1,V2));
end

c=299792458;
[an,bn,RCSTheta1D,ETheta1D] = mieHKURCS(rad,c/lambda,n_imm^2,1,nsphere^2,1,40,theta1D);
ETheta = reshape(ETheta1D, [size(X,1) size(X,2)]);

R = sqrt(X.^2 + Y.^2 + z(end)^2);
E_sca_Mie = ETheta./( exp(1i*k*z(end))/z(end) ).*( exp(1i*k*R)./R );
E_sca_Mie = E_sca_Mie./E_sca_Mie(end/2,end/2);


NA_array = [0.2790, 0.4776, 0.6584];
N_NA = size(NA_array, 2);

%x = x./lambda;
%y = y./lambda;
%rad = rad/lambda;

for jj = 1:size(defocus,2)

    Z = -L(3)/2+defocus(jj);

for ii = 1:N_NA

NA =NA_array(ii);

E_MLB_sca = BPM(E_sca_MLB, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MLR_sca = BPM(E_sca_MLR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR_sca = BPM(E_sca_MSR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSL_sca = BPM(E_sca_MSL, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_Mie_sca = BPM(E_sca_Mie, Z, NA, k, 2*pi*fxx, 2*pi*fyy);

E_MLB_tot = BPM(E_tot_MLB, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MLR_tot = BPM(E_tot_MLR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR_tot = BPM(E_tot_MSR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSL_tot = BPM(E_tot_MSL, Z, NA, k, 2*pi*fxx, 2*pi*fyy);

%E_fft_last = fftshift(fft2(E_tot_MLB(:,:,end)));
%E_fft_mid  = fftshift(fft2(E_MLB_tot));
%figure
%subplot(2, 2, 1)
%imagesc(fx,fy,abs(E_fft_last))
%axis square
%colorbar
%title('E_{last}')

%subplot(2, 2, 2)
%imagesc(fx,fy,abs(E_fft_mid))
%axis square
%colorbar
%title('E_{mid}')

%subplot(2, 2, 3)
%imagesc(fx,fy,abs(E_tot_MLB(:,:,end)))
%axis square
%colorbar
%title('E_{last}')

%subplot(2, 2, 4)
%imagesc(fx,fy,abs(E_MLB_tot))
%axis square
%colorbar
%title('E_{mid}')

cmax_sca_E = max(max(abs(E_MLB_sca)));
cmax_sca_I = max(max(abs(E_MLB_sca).^2));
cmax_tot_E = max(max(abs(E_MLB_tot)));
cmax_tot_I = max(max(abs(E_MLB_tot).^2));
cmin_tot_E = min(min(abs(E_MLB_tot)));
cmin_tot_I = min(min(abs(E_MLB_tot).^2));

figure
set(gcf, 'Position', get(0, 'Screensize'));
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


angle_MLB = unwrap(angle(E_MLB_sca(end/2,:)));
angle_MLR = unwrap(angle(E_MLR_sca(end/2,:)));
angle_MSR = unwrap(angle(E_MSR_sca(end/2,:)));
angle_Mie = unwrap(angle(E_Mie_sca(end/2,:)));
% remornalize
angle_MLB = angle_MLB - angle_MLB(end/2);
angle_MLR = angle_MLR - angle_MLR(end/2);
angle_MSR = angle_MSR - angle_MSR(end/2);
angle_Mie = angle_Mie - angle_Mie(end/2);

subplot(5, 6, [28 29 30])
hold on
plot(x, angle_MLB,'LineWidth',4);
plot(x, angle_MLR,'LineStyle',"--",'LineWidth',4);
plot(x, angle_MSR,'LineStyle',":",'LineWidth',4);
plot(x, angle_Mie,'LineStyle',"-.",'LineWidth',4);
legend('MLB','MLR','MSR','Mie')
xlim([-2*rad 2*rad])

f_title = sprintf('d = %dum, NA_{in} = %1.4f, NA_{obj} = %1.4f, defocus = %dum', 2*rad, NA_in, NA_array(ii), defocus(jj));
sgtitle(f_title)
%saveas(gcf,'Figures/NA0p0.png')
saveas(gcf,[pwd sprintf('/Figures/d%dNA0p%dNAimag%ddefocus%d.jpg',round(rad*2),round(10*NA_in),round(NA_array(ii)),round(defocus(jj)))]);

end

end






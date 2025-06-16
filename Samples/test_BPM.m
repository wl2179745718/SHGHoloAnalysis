clear variables; close all; clc; addpath(genpath('../Functions'));

ps          = .105/2;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.5;                  % central wavelength (microns)  
n_imm       = 1;                % refractive index of immersion media
nsphere=1.02;%1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.8;
NA_in = round(NA_in/deltaf(1)/lambda)*deltaf(1)*lambda;
phi       = asin(NA_in);              

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
rad = 2*lambda;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=1;

U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

ord = 1;

E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017

V_inc = 0*V;
E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,'Vol');
E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
U_inp_end_MSR = E_MSR_3d_inc(:,:,end);
U_inp_end_MLR = E_MLR_3d_inc(:,:,end);
U_inp_end_MLB = E_MLB_3d_inc(:,:,end);

figure
subplot(2,3,1)
imagesc(x,z,abs(U_inp_end_MLB))
axis square
colorbar
subplot(2,3,4)
imagesc(x,z,angle(U_inp_end_MLB))
axis square
colorbar
subplot(2,3,2)
imagesc(x,z,abs(U_inp_end_MLR))
axis square
colorbar
subplot(2,3,5)
imagesc(x,z,angle(U_inp_end_MLR))
axis square
colorbar
subplot(2,3,3)
imagesc(x,z,abs(U_inp_end_MSR))
axis square
colorbar
subplot(2,3,6)
imagesc(x,z,angle(U_inp_end_MSR))
axis square
colorbar

figure
subplot(2,3,1)
imagesc(x,z,abs(E_MLB_3d(:,:,end)))
axis square
colorbar
subplot(2,3,4)
imagesc(x,z,angle(E_MLB_3d(:,:,end)))
axis square
colorbar
subplot(2,3,2)
imagesc(x,z,abs(E_MLR_3d(:,:,end)))
axis square
colorbar
subplot(2,3,5)
imagesc(x,z,angle(E_MLR_3d(:,:,end)))
axis square
colorbar
subplot(2,3,3)
imagesc(x,z,abs(E_MSR_3d(:,:,end)))
axis square
colorbar
subplot(2,3,6)
imagesc(x,z,angle(E_MSR_3d(:,:,end)))
axis square
colorbar

figure
subplot(2,3,1)
imagesc(x,z,abs(E_MLB_3d(:,:,end)-U_inp_end_MLB))
axis square
colorbar
subplot(2,3,4)
imagesc(x,z,angle(E_MLB_3d(:,:,end)-U_inp_end_MLB))
axis square
colorbar
subplot(2,3,2)
imagesc(x,z,abs(E_MLR_3d(:,:,end)-U_inp_end_MLR))
axis square
colorbar
subplot(2,3,5)
imagesc(x,z,angle(E_MLR_3d(:,:,end)-U_inp_end_MLR))
axis square
colorbar
subplot(2,3,3)
imagesc(x,z,abs(E_MSR_3d(:,:,end)-U_inp_end_MSR))
axis square
colorbar
subplot(2,3,6)
imagesc(x,z,angle(E_MSR_3d(:,:,end)-U_inp_end_MSR))
axis square
colorbar

E_MSR_xz = squeeze(E_MSR_3d(end/2,:,:));
E_MLR_xz = squeeze(E_MLR_3d(end/2,:,:));
E_MLB_xz = squeeze(E_MLB_3d(end/2,:,:));
figure
subplot(2,3,1)
imagesc(x,z,abs(E_MLB_xz))
axis square
colorbar
subplot(2,3,4)
imagesc(x,z,angle(E_MLB_xz))
axis square
colorbar
subplot(2,3,2)
imagesc(x,z,abs(E_MLR_xz))
axis square
colorbar
subplot(2,3,5)
imagesc(x,z,angle(E_MLR_xz))
axis square
colorbar
subplot(2,3,3)
imagesc(x,z,abs(E_MSR_xz))
axis square
colorbar
subplot(2,3,6)
imagesc(x,z,angle(E_MSR_xz))
axis square
colorbar

[E_tot_MSR,E_sca_MSR] = tot2sca(E_MSR_3d,U_inp_end_MSR);
[E_tot_MLR,E_sca_MLR] = tot2sca(E_MLR_3d,U_inp_end_MLR);
[E_tot_MLB,E_sca_MLB] = tot2sca(E_MLB_3d,U_inp_end_MLB);

figure
cmax = max(max(abs(E_tot_MSR)));
subplot(2,3,1)
imagesc(x,y,abs(E_tot_MSR))
axis square
clim([0 cmax]);
colorbar
title('MSR |E_{tot}|')

subplot(2,3,4)
imagesc(x,y,abs(E_sca_MSR))
clim([0 cmax]);
axis square
colorbar
title('MSR |E_{sca}|')

subplot(2,3,2)
imagesc(x,y,abs(E_tot_MLR))
axis square
clim([0 cmax]);
colorbar
title('MLR |E_{tot}|')

subplot(2,3,5)
imagesc(x,y,abs(E_sca_MLR))
clim([0 cmax]);
axis square
colorbar
title('MLR |E_{sca}|')

subplot(2,3,3)
imagesc(x,y,abs(E_tot_MLB))
axis square
clim([0 cmax]);
colorbar
title('MLB |E_{tot}|')

subplot(2,3,6)
imagesc(x,y,abs(E_sca_MLB))
clim([0 cmax]);
axis square
colorbar
title('MLB |E_{sca}|')

E_sca_MSR = E_sca_MSR./E_sca_MSR(end/2,end/2);
E_sca_MLR = E_sca_MLR./E_sca_MLR(end/2,end/2);
E_sca_MLB = E_sca_MLB./E_sca_MLB(end/2,end/2);

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


Z = -L(3)/2;

NA_array = [0.2, 0.4, 0.8];
N_NA = size(NA_array, 2);

x = x./lambda;
y = y./lambda;
rad = rad/lambda;

for ii = 1:N_NA

NA =NA_array(ii);

E_MLB = BPM(E_sca_MLB, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MLR = BPM(E_sca_MLR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR = BPM(E_sca_MSR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_Mie = BPM(E_sca_Mie, Z, NA, k, 2*pi*fxx, 2*pi*fyy);

figure(100)
set(gcf, 'Position', get(0, 'Screensize'));
colormap(gray(256));
subplot(5, 2*N_NA, 2*(ii-1)+1)
imagesc(x,y,abs(E_MLB))
axis square
colorbar
title('MLB |E|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 2*ii)
imagesc(x,y,angle(E_MLB))
axis square
colorbar
title('MLB \angle E')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 2*N_NA+2*(ii-1)+1)
imagesc(x,y,abs(E_MLR))
axis square
colorbar
title('MLR |E|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 2*N_NA+2*ii)
imagesc(x,y,angle(E_MLR))
axis square
colorbar
title('MLR \angle E')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 4*N_NA+2*(ii-1)+1)
imagesc(x,y,abs(E_MSR))
axis square
colorbar
title('MSR |E|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 4*N_NA+2*ii)
imagesc(x,y,angle(E_MSR))
axis square
colorbar
title('MSR \angle E')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 6*N_NA+2*(ii-1)+1)
imagesc(x,y,abs(E_Mie))
axis square
colorbar
title('Mie |E|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(5, 2*N_NA, 6*N_NA+2*ii)
imagesc(x,y,angle(E_Mie))
axis square
colorbar
title('Mie \angle E')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

angle_MLB = unwrap(angle(E_MLB(end/2,:)));
angle_MLR = unwrap(angle(E_MLR(end/2,:)));
angle_MSR = unwrap(angle(E_MSR(end/2,:)));
angle_Mie = unwrap(angle(E_Mie(end/2,:)));
% remornalize
angle_MLB = angle_MLB - angle_MLB(end/2);
angle_MLR = angle_MLR - angle_MLR(end/2);
angle_MSR = angle_MSR - angle_MSR(end/2);
angle_Mie = angle_Mie - angle_Mie(end/2);

subplot(5, 2*N_NA, [8*N_NA+2*(ii-1)+1 8*N_NA+2*ii])
hold on
plot(x, angle_MLB,'LineWidth',4);
plot(x, angle_MLR,'LineStyle',"--",'LineWidth',4);
plot(x, angle_MSR,'LineStyle',":",'LineWidth',4);
plot(x, angle_Mie,'LineStyle',"-.",'LineWidth',4);
legend('MLB','MLR','MSR','Mie')
xlim([-2*rad 2*rad])

end

figure(100)
f_title = sprintf('NA_{in} = %1.2f, NA_{obj} = %1.2f, %1.2f, %1.2f', NA_in, NA_array(1), NA_array(2),NA_array(3));
sgtitle(f_title)
%saveas(gcf,'Figures/NA0p0.png')
saveas(gcf,[pwd sprintf('/Figures/NA0p%d.jpg',round(10*NA_in))]);


figure('units','normalized','outerposition',[0 0 1 1])
%cmax = max(max(abs(E_sca_MSR)));
cmax2 = max(max(abs(E_MSR)));

subplot(4,6,1)
imagesc(x,y,abs(E_sca_MLB))
axis square
clim([0 cmax]);
colormap(gray(256));
colorbar
title('MLB |E_{sca}|')

subplot(4,6,2)
imagesc(x,y,angle(E_sca_MLB))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,6,3)
imagesc(x,y,abs(E_MLB))
axis square
clim([0 cmax2]);
colorbar
title('MLB |E_{sca}|')

subplot(4,6,4)
imagesc(x,y,angle(E_MLB))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,6,5)
imagesc(x,y,abs(E_MLB))
axis square
clim([0 cmax2]);
colorbar
title('MLB |E_{sca}|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,6)
imagesc(x,y,angle(E_MLB))
axis square
colorbar
title('\angle E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,7)
imagesc(x,y,abs(E_sca_MSR))
axis square
clim([0 cmax]);
colorbar
title('MSR |E_{sca}|')

subplot(4,6,8)
imagesc(x,y,angle(E_sca_MSR))
axis square
colorbar
title('MSR \angle E_{sca}')

subplot(4,6,9)
imagesc(x,y,abs(E_MSR))
axis square
clim([0 cmax2]);
colorbar
title('MSR |E_{sca}|')

subplot(4,6,10)
imagesc(x,y,angle(E_MSR))
axis square
colorbar
title('MSR \angle E_{sca}')

subplot(4,6,11)
imagesc(x,y,abs(E_MSR))
axis square
clim([0 cmax2]);
colorbar
title('MSR |E_{sca}|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,12)
imagesc(x,y,angle(E_MSR))
axis square
colorbar
title('MSR \angle E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,13)
imagesc(x,y,abs(E_sca_MLR))
axis square
clim([0 cmax]);
colorbar
title('MLR |E_{sca}|')

subplot(4,6,14)
imagesc(x,y,angle(E_sca_MLR))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,6,15)
imagesc(x,y,abs(E_MLR))
axis square
clim([0 cmax2]);
colorbar
title('MLR |E_{sca}|')

subplot(4,6,16)
imagesc(x,y,angle(E_MLR))
axis square
colorbar
title('\angle E_{sca}')

subplot(4,6,17)
imagesc(x,y,abs(E_MLR))
axis square
clim([0 cmax2]);
colorbar
title('MLR |E_{sca}|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,18)
imagesc(x,y,angle(E_MLR))
axis square
colorbar
title('\angle E_{sca}')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,19)
imagesc(x,y,abs(E_sca_Mie))
clim([0 cmax]);
axis square
colorbar
title('ground truth |Mie|')

subplot(4,6,20)
imagesc(x,y,angle(E_sca_Mie))
axis square
colorbar
title('ground truth \angle Mie')

subplot(4,6,21)
imagesc(x,y,abs(E_Mie))
clim([0 cmax2]);
axis square
colorbar
title('ground truth |Mie|')

subplot(4,6,22)
imagesc(x,y,angle(E_Mie))
axis square
colorbar
title('ground truth \angle Mie')

subplot(4,6,23)
imagesc(x,y,abs(E_Mie))
clim([0 cmax2]);
axis square
colorbar
title('ground truth |Mie|')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

subplot(4,6,24)
imagesc(x,y,angle(E_Mie))
axis square
colorbar
title('ground truth \angle Mie')
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

angle_MLB = unwrap(angle(E_MLB(end/2,:)));
angle_MLR = unwrap(angle(E_MLR(end/2,:)));
angle_MSR = unwrap(angle(E_MSR(end/2,:)));
angle_Mie = unwrap(angle(E_Mie(end/2,:)));
% remornalize
angle_MLB = angle_MLB - angle_MLB(end/2);
angle_MLR = angle_MLR - angle_MLR(end/2);
angle_MSR = angle_MSR - angle_MSR(end/2);
angle_Mie = angle_Mie - angle_Mie(end/2);

figure('units','normalized','outerposition',[0 0 1 1])
hold on
plot(x, angle_MLB,'LineWidth',4);
plot(x, angle_MLR,'LineStyle',"--",'LineWidth',4);
plot(x, angle_MSR,'LineStyle',":",'LineWidth',4);
plot(x, angle_Mie,'LineStyle',"-.",'LineWidth',4);
legend('MLB','MLR','MSR','Mie')
xlim([-2*rad 2*rad])

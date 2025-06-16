clear variables; close all; clc; addpath(genpath('../Functions'));

ps          = .105;                 % pixel size (x,y,z) in object space (microns)
lambda      = 1;                  % central wavelength (microns)
NA_in = 0.2;   
n_imm       = 1;                % refractive index of immersion media
nsphere=1.00;%1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
deltaf = 1./L;
Lf = 1./delta;
[fx,fy,fz] = L2xyz(Lf,deltaf);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

NA_in = round(NA_in/deltaf(1))*deltaf(1);
phi       = asin(NA_in);             

dGk = 0;
rad = 2*lambda;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=1;

U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

E_MLB_3d=MultiLayerBornv3(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');

th = 0:pi/50:2*pi;
xunit = k/2/pi * cos(th);
yunit = k/2/pi * sin(th);

holo_video = VideoWriter('incident','MPEG-4'); %create the video object
holo_video.FrameRate = 5;
open(holo_video);

%x = x./lambda;
%y = y./lambda;
%z = z./lambda;
for ii = 1%:10:N(3)

E_xy = E_MLB_3d(:,:,ii);
E_kxy = fftshift(fft2(E_xy));

figure
set(gcf, 'Position', get(0, 'Screensize'));
colormap(gray(256));
subplot(2,2,1)
imagesc(x,y,abs(E_xy))
clim([0 1.2])
axis square
colorbar
subplot(2,2,2)
imagesc(x,y,angle(E_xy))
clim([-pi pi])
axis square
colorbar
subplot(2,2,3)
imagesc(fx,fy,abs(E_kxy))
hold on
plot(xunit, yunit);
xlim([-k/2/pi*1.5 k/2/pi*1.5])
ylim([-k/2/pi*1.5 k/2/pi*1.5])
clim([0 2.4e5])
axis square
colorbar
subplot(2,2,4)
imagesc(fx,fy,angle(E_kxy))
hold on
plot(xunit, yunit);
xlim([-k/2/pi*1.5 k/2/pi*1.5])
ylim([-k/2/pi*1.5 k/2/pi*1.5])
clim([-pi pi])
axis square
colorbar

f_title = sprintf('z = %1.2f', z(ii));
sgtitle(f_title)

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20.8 10.8]);
saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file

end

close(holo_video);


E_MLB_xz = squeeze(E_MLB_3d(end/2,:,:));
figure
subplot(1,2,1)
imagesc(x,y,abs(E_MLB_xz))
axis square
colorbar
subplot(1,2,2)
imagesc(x,y,angle(E_MLB_xz))
axis square
colorbar

clear variables; close all; clc; addpath(genpath('../Functions')); addpath(genpath('../Data&Tutorial'));
init_unlocbox();

data_origin      = '/Users/langwang/Documents/GitHub/SHGHoloAnalysis/Data&Tutorial/CZI scattering data-selected/LowNA';
save_recon       = true;

sample = load([data_origin '/700_sample.mat'],'acq');
bkgd = load([data_origin '/700_background.mat'],'acq');
disp('finished loading measurements');

load([data_origin '/700_k_illumination_calibrated.mat']);
disp('finished loading computationally-corrected NA and illumination-angle parameters');

dwnsmple_fac        = 1; % manual correction in case raw measurement is oversampled (downsample_fac=1 means no downsampling)
totFOV_amp_acqs     = sample.acq ./ bkgd.acq;
totFOV_amp_acqs     = imresize(totFOV_amp_acqs,dwnsmple_fac);
size(totFOV_amp_acqs)

disp('finished loading data and data_parameters');

ps          = .105/2;                 % pixel size (x,y,z) in object space (microns)
lambda      = 0.7;                  % central wavelength (microns)  
n_imm       = 1.5678;                % refractive index of immersion media
nsphere=1.5717;%1.2;
n=[nsphere,n_imm];
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
rad = 7/2;
NA = 0.2790 ;
Z = -L(3)/2;

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);
ord = 1;

Eps=n_imm^2/lambda^2*0.25;

holo_video = VideoWriter('exp_low_4','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);

test = asin(sqrt(kx_corr.^2+ky_corr.^2)./(k/2/pi))./pi.*180;

for ii = 1:size(totFOV_amp_acqs,3)
ii
fx_in           = kx_corr(ii);                               
fy_in           = ky_corr(ii);  

fx_in = round(fx_in/deltaf(1))*deltaf(1);
fy_in = round(fy_in/deltaf(2))*deltaf(2);
%fx_in_interp = fix(fx_in/dfx)*dfx;
%fy_in_interp = fix(fy_in/dfx)*dfx;
U_inp         = exp(1i * 2 * pi * (fx_in * X + fy_in * Y));

E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MSR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
E_MSL_3d=MultiSlice(U_inp,RI,lambda,ps,fxx,fyy,n_imm);
E_tot_MLB = E_MLB_3d(:,:,end);
E_tot_MSR = E_MSR_3d(:,:,end);
E_tot_MSL = E_MSL_3d(:,:,end);

E_MLB_tot = BPM(E_tot_MLB, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSR_tot = BPM(E_tot_MSR, Z, NA, k, 2*pi*fxx, 2*pi*fyy);
E_MSL_tot = BPM(E_tot_MSL, Z, NA, k, 2*pi*fxx, 2*pi*fyy);

tiledlayout(1,4,'TileSpacing','none');
%set(gcf, 'Position', get(0, 'Screensize'));
nexttile
imagesc(totFOV_amp_acqs(674:734,598:658,ii));  % Low
%imagesc(totFOV_amp_acqs(664:724,604:664,ii));  % Med
%imagesc(totFOV_amp_acqs(658:718,630:690,ii));  % High
axis equal; colormap gray; axis tight; clim(1+[-1.0,1.4]); set(gca,'YDir','normal'); axis off

nexttile
imagesc(x,y,abs(E_MLB_tot).^2)
axis equal; colormap gray; axis tight; set(gca,'YDir','normal'); axis off
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

nexttile
imagesc(x,y,abs(E_MSR_tot).^2)
axis equal; colormap gray; axis tight; set(gca,'YDir','normal'); axis off
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])

nexttile
imagesc(x,y,abs(E_MSL_tot).^2)
axis equal; colormap gray; axis tight; set(gca,'YDir','normal'); axis off
xlim([-2*rad 2*rad])
ylim([-2*rad 2*rad])


saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file


end

close(holo_video);
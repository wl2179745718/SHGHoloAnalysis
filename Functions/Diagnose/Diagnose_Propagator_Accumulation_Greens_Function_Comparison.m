function [] = Diagnose_Propagator_Accumulation_Greens_Function_Comparison(X, Y, z_positive, FX, FY, fx, fy, dr, k0, lambda)

%G2 = -exp(1i*k0*sqrt((dr(3))^2+ X.^2 + Y.^2))./sqrt((dr(3))^2+ X.^2 + Y.^2); 
%G2F = myfft(G2);
%G2F( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;

n_imm=1;
dGk = -0.001;
Eps = 0;
%G2k = G_kx_ky(FX, FY,n_imm,lambda,dr(3),dGk,Eps);
%G2k( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;


holo_video = VideoWriter('Greens function comparison','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);


for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation to compare 2 Greens function: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2)/4/pi; 
G1F = myfft(G1);
G1F( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;

G2F = G_kx_ky(FX, FY,n_imm,lambda,z_positive(ii),dGk,Eps)./dr(1)./dr(2);
G2F( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;

errorG(ii) = norm(G1F - G2F,2)/norm(G1F,2);

figure
subplot(1,3,1)
imagesc(lambda*fx,lambda*fy,abs(G1F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('clipped F\{ G_r \}')

subplot(1,3,2)
imagesc(lambda*fx,lambda*fy,abs(G2F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('clipped  G_k')

subplot(1,3,3)
imagesc(lambda*fx,lambda*fy,abs(G1F - G2F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('difference')

set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file


end

figure('Name','Diagnose: compare Greens functions in r and k space')
subplot(2,1,1)
plot(z_positive/lambda,errorG)
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Greens function comparison.png");
imshow(A,map)

end
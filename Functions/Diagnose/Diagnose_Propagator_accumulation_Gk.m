function [] = Diagnose_Propagator_accumulation_Gk(z_positive, FX, FY, fx, fy, dr, Pdz, lambda)

%G2 = -exp(1i*k0*sqrt((dr(3))^2+ X.^2 + Y.^2))./sqrt((dr(3))^2+ X.^2 + Y.^2); 
n_imm=1;
dGk = -0.001;
Eps = 0;
G2F = G_kx_ky(FX, FY,n_imm,lambda,dr(3),dGk,Eps);

%G2 = myifft(G2F);


holo_video = VideoWriter('clipped Greens function in k space','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);


for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation with clipped Greens function in k space: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1F = G_kx_ky(FX, FY,n_imm,lambda,z_positive(ii),dGk,Eps);

errorG(ii) = norm(G1F - G2F,2)/norm(G1F,2);

figure
subplot(1,3,1)
imagesc(lambda*fx,lambda*fy,abs(G1F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,2)
imagesc(lambda*fx,lambda*fy,abs(G2F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,3)
imagesc(lambda*fx,lambda*fy,abs(G1F - G2F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file


G2F = Pdz.*(G2F); 


end

figure('Name','Diagnose: propagator error accumulation with clipped Greens function in k space')
subplot(2,1,1)
plot(z_positive/lambda,errorG)
text(z_positive(end)/lambda/2, 0.9*max(errorG), ['dz = ', num2str(dr(3)/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Clipped Greens Function in k space.png");
imshow(A,map)

end
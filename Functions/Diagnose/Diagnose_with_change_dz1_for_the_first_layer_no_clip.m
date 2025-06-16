function [] = Diagnose_with_change_dz1_for_the_first_layer_no_clip(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda, dz1)

G2 = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2); 
G2F = myfft(G2);
barMax = max(max( log(abs(G2F)) ));
barMin = log(abs(G2F(round(end/2),round(end/2))))-3;
G3_0 = G2F;

holo_video = VideoWriter('unclipped Greens function','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);

z_positive = z_positive - dr(3) + dz1;

for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation with clipped Greens function: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2); 
G1F = myfft(G1);

errorG(ii) = norm(G1 - G2,2)/norm(G1,2);

figure
subplot(1,3,1)
imagesc(lambda*fx,lambda*fy,log(abs(G1F)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
clim([barMin,barMax])
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('log(|G1_k|)')

subplot(1,3,2)
imagesc(lambda*fx,lambda*fy,log(abs(myfft(G2))))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
clim([barMin,barMax])
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('log(|G2_k|)')

subplot(1,3,3)
imagesc(lambda*fx,lambda*fy,log(abs(myfft(G1 - G2))))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
clim([barMin,barMax])
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])
title('log(|G1_k-G2_k|)')

global_title=['z= ', num2str(z_positive(ii))];
sgtitle(global_title)
set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file


G2 = myifft(Pdz.*(myfft(G2))); 
G3 = myifft(Propagator(lambda,FX,FY,z_positive(ii)).*(G3_0));

errorG2G3(ii) = norm(G2-G3,2)/norm(G2,2);

end

figure('Name','Diagnose: propagator error accumulation with clipped Greens function')
subplot(2,1,1)
plot(z_positive/lambda,errorG)
text(z_positive(end)/lambda/2, 0.2*max(errorG)+0.8*min(errorG), ['dz1 = ', num2str(dz1/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose with changing dz for the first layer no clip.png");
imshow(A,map)

figure('Name','Diagnose: error between G2 and G3')
plot(z_positive,errorG2G3)

end
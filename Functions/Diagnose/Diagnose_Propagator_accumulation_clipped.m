function [] = Diagnose_Propagator_accumulation_clipped(X, Y, z_positive, FX, FY, fx, fy, dr, Pdz, k0, lambda)

G2 = -exp(1i*k0*sqrt((dr(3))^2+ X.^2 + Y.^2))./sqrt((dr(3))^2+ X.^2 + Y.^2); 
G2F = myfft(G2);
%G2F( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;
G2 = myifft(G2F);
G3_0 = G2F;

holo_video = VideoWriter('clipped Greens function','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);


for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation with clipped Greens function: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2); 
G1F = myfft(G1);
%G1F( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;
G1 = myifft(G1F);

errorG(ii) = norm(G1 - G2,2)/norm(G1,2);

figure
subplot(1,3,1)
imagesc(lambda*fx,lambda*fy,abs(G1F))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,2)
imagesc(lambda*fx,lambda*fy,abs(myfft(G2)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,3)
imagesc(lambda*fx,lambda*fy,abs(myfft(G1 - G2)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

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
text(z_positive(end)/lambda/2, 0.9*max(errorG), ['dz = ', num2str(dr(3)/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Clipped Greens Function in k space2.png");
imshow(A,map)

figure('Name','Diagnose: error between G2 and G3')
plot(z_positive,errorG2G3)

end
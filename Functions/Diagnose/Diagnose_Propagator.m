function [] = Diagnose_Propagator(X, Y, z, dr, x0_index,y0_index, z0_index, fx, fy, Pdz, k0, lambda)

src = zeros(size(X)); 

src(x0_index,y0_index)=1; 

z_positive = z(z0_index+1:end);


holo_video = VideoWriter('propagator error','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);


for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii)+dr(3))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii)+dr(3))^2+ X.^2 + Y.^2); 

G2 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2); 

G3 = myifft(Pdz.*(myfft(G2))); 

errorG(ii) = norm(G1 - G3,2)/norm(G1,2);


figure
subplot(1,3,1)
imagesc(lambda*fx,lambda*fy,abs(myfft(G1)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,2)
imagesc(lambda*fx,lambda*fy,abs(myfft(G3)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(1,3,3)
imagesc(lambda*fx,lambda*fy,abs(myfft(G1 - G3)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file

end

figure('Name','Diagnose: error in the propagator')
subplot(2,1,1)
plot(z_positive/lambda,errorG)
text(z_positive(end)/lambda/2, 0.9*max(errorG), ['dz = ', num2str(dr(3)/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose Propagator.png");
imshow(A,map)



end
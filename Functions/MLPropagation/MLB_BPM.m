function [Us_xy, Us_xz, Us_yz]=MLB_BPM(shape,U_in,AG,propdz, uin_xz, uin_yz, x0_index, y0_index)%, test_V, test_Us, test_AG

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);

Us_xz = zeros(size(uin_xz));
Us_yz = zeros(size(uin_yz));

%V_tst = zeros(size(U_in));

holo_video = VideoWriter('back-propagated image along z','MPEG-4'); %create the video object
holo_video.FrameRate = 10;
open(holo_video);

U = U_in; % Initial Field
for i=1:Box_N(3)
    S=U;
    U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    start=[1 1 i]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    Us=myifft((myfft(S.*Vn)).*(AG))*Box_delta(1)*Box_delta(2)*Box_delta(3); % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    %if sum(Vn,'all')~=0
    %    test_V = Vn;
    %    test_Us = Us;
    %    test_AG = AG;
    %end
    %U_t(i) = S(x0_index, y0_index);

    U=U+Us;

    %V_tst = V_tst+Vn*Box_delta(1)*Box_delta(2)*Box_delta(3);

    Us_yz(i,:) = U(:,x0_index).' - uin_yz(i,:);
    Us_xz(i,:) = U(y0_index,:) - uin_xz(i,:);

    uin_xy = exp(1i*2*pi*fxin*X)*exp(1i*2*pi*(z(i)-z(1)+dr(3))*cos(thetain)/lambda);

    Us_xy = U - uin_xy;

    PZ = Propagator(lambda,FX,FY, -z(i));
    E_MLB_image_k = PZ.*(myfft(Us_xy));
    E_MLB_image=myifft(E_MLB_image_k);

    figure
subplot(2,2,1)
imagesc(lambda*fx,lambda*fy,abs(E_MLB_image_k))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(2,2,3)
imagesc(lambda*fx,lambda*fy,wrapToPi(angle(E_MLB_image_k)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(2,2,2)
imagesc(x,y,abs(E_MLB_image))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

subplot(2,2,4)
imagesc(x,y,wrapToPi(angle(E_MLB_image)))
axis equal; colormap gray; axis tight; colorbar%clim(1+[-1.0,1.4]);
set(gca,'YDir','normal'); %axis off
xlim([-2 2])
ylim([-2 2])

sgtitle('z = ')
set(gcf, 'Position', get(0, 'Screensize'));

saveas(gcf,'Figure.png');
thisimage = imread('Figure.png');
writeVideo(holo_video,thisimage); %write the image to file


end

end
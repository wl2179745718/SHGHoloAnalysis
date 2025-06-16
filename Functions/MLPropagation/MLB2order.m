function [Us_xy, Us_xz, Us_yz] = MLB2order(shape,U_in,AG,propdz,prop2dz, uin_xy, uin_xz, uin_yz, x0_index, y0_index)

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);

U = U_in;

Us_xz = zeros(size(uin_xz));
Us_yz = zeros(size(uin_yz));

for i=1:2:Box_N(3)
    S=U;
    U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    start=[1 1 i]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    Us=myifft((myfft(S.*Vn)).*(AG))*Box_delta(1)*Box_delta(2)*Box_delta(3); % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
    Us_yz(i,:) = U(:,x0_index).' - uin_yz(i,:);
    Us_xz(i,:) = U(y0_index,:) - uin_xz(i,:);

    % The second half of the layer
    start=[1 1 i+1]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    Us=myifft((myfft(U.*Vn)).*(AG))*2*Box_delta(1)*Box_delta(2)*Box_delta(3);
    U = Us + myifft(prop2dz.*(myfft(S)));
    Us_yz(i+1,:) = U(:,x0_index).' - uin_yz(i+1,:);
    Us_xz(i+1,:) = U(y0_index,:) - uin_xz(i+1,:);

end

Us_xy = U - uin_xy;

end
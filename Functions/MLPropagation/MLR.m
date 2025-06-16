function [Us_xy, Us_xz, Us_yz]=MLR(shape,U_in,AG,propdz, uin_xy, uin_xz, uin_yz, x0_index, y0_index)


% This function computes the propagation of an initial field U_in through a
% scattering potential V using a first Rytov approximation for each slice of V. 

% Inputs: fxx,fyy: Spatial frequency grids (ifftshifted)
%         lambda: wavelength in microns
%         n_imm: index of refraction of background material
%         dz:    Thickness of layer in microns
%         V:     Scattering potential of object/material 3d matrix array
%         U_in:  input field
%         eps:   absorption term
%         opt:   options either 'out' or 'Vol'. 'out' just outputs the
%         final output field. 'Vol' outs a 3d array containing the field as it propagates through each
%         slice

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);

%AG = G;
%propdz = P;
U=(U_in); % Initial Field

Us_xz = zeros(size(uin_xz));
Us_yz = zeros(size(uin_yz));

%S=zeros(size(V));
for i=1:Box_N(3)
    S=U;
    U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
    start=[1 1 i]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    Us=myifft((myfft(S.*Vn)).*(AG))*Box_delta(1)*Box_delta(2)*Box_delta(3); % Scattered field of nth layer
    %Us2=ifft2((fft2(Us.*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz;
    %U=U+Us;
    %U=U.*exp((Us./U)+(Us2./U)-.5*(Us2./U).^2);
    B1=Us./U;
    B1(abs(B1)>1)=0; % Avoiding indeterminants/large numbers (from zeros in field) assumes phase is small
%     B2=Us2./U;
%     B2(abs(B2)>1)=0; % Added these to avoid indeterminants/large numbers
    %U=U.*exp((B1)+(B2)-.5*((B2)).^2);
    U=U.*exp((B1));
    Us_yz(i,:) = U(:,x0_index).' - uin_yz(i,:);
    Us_xz(i,:) = U(y0_index,:) - uin_xz(i,:);
end

Us_xy = U - uin_xy;

end
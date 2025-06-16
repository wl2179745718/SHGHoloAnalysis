function U=MultiLayerBornv3(shape,dz,inci_flag,U_in,AG,propdz)

% This function computes the propagation of an initial field U_in through a
% scattering potential V using a first Born approximation for each slice of V. 

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

%prop_crop= (fxx.^2 + fyy.^2 > (n_imm/lambda)^2==0); Elimination of
%evanescent waves
filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [sphere_N(1),  sphere_N(2)];
Nz = sphere_N(3);
%AG = G;
%propdz = P;
U = U_in; % Initial Field
for i=1:Nz-1
    S=U;
    U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    if inci_flag == 1
        Vn = zeros(chunk_size);
    else
        start=[1 1 i]; % indicates which layer to read from the data file
        count=[chunk_size 1]; % Chunk size
        Vn = h5read('../Medium/sphere.h5','/../Medium/sphere',start,count);
    end
    Us=(ifft2((fft2(S.*Vn)).*(AG))*dz); % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
end

end
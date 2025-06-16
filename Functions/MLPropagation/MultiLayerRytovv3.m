function U=MultiLayerRytovv3(shape,dz,inci_flag,U_in,AG,propdz)


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
chunk_size = [sphere_N(1),  sphere_N(2)];
Nz = sphere_N(3);

%AG = G;
%propdz = P;
U=(U_in); % Initial Field

%S=zeros(size(V));
for i=1:Nz
    S=U;
    U=ifft2(propdz.*(fft2(U))); % Prop dz to next layer
    if inci_flag == 1
        Vn = zeros(chunk_size);
    else
        start=[1 1 i]; % indicates which layer to read from the data file
        count=[chunk_size 1]; % Chunk size
        Vn = h5read('../Medium/sphere.h5','/../Medium/sphere',start,count);
    end
    Us=ifft2((fft2(S.*Vn)).*(AG))*dz; % Scattered field of nth layer
    %Us2=ifft2((fft2(Us.*V(:,:,i))).*(AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,dz)))*dz;
    %U=U+Us;
    %U=U.*exp((Us./U)+(Us2./U)-.5*(Us2./U).^2);
    B1=Us./U;
    B1(abs(B1)>1)=0; % Avoiding indeterminants/large numbers (from zeros in field) assumes phase is small
%     B2=Us2./U;
%     B2(abs(B2)>1)=0; % Added these to avoid indeterminants/large numbers
    %U=U.*exp((B1)+(B2)-.5*((B2)).^2);
    U=U.*exp((B1));
end

end
function evol=MultiLayerBornv2(fxx,fyy,lambda,n_imm,dz,V,U_in,Eps,dGk,opt)

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
bound = dGk * ( max(max(sqrt(fxx.^2 + fyy.^2))) - n_imm/lambda );

SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field

prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2));%+1i*eps %Add small absorption term to avoid indeterminates in the angular greens function


%prop=@(z) exp(prop_phs*z);
%Mask=(n_imm/lambda)^2>1.01*((fxx.^2+fyy.^2));

AG= ((-1i.*exp(prop_phs.*dz)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps)))); % Angular Greens function
AG(fxx.^2 + fyy.^2 > (n_imm/lambda + bound)^2) = 0;

%AG(isnan(AG)==1)=0;
propdz=exp(dz*prop_phs);%prop(dz);
U=(U_in); % Initial Field
switch opt
    case 'out'
        k0=(2*pi)/lambda;
for i=1:V(2)
    S=U;
    U=ifft2(propdz.*(fft2(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    if V(1) == 0
        Vn = zeros(size(fxx));
    else
        start=[1 1 i]; % indicates which layer to read from the data file
        count=[size(fxx) 1]; % Chunk size
        RI = h5read('sphere.h5','/sphere',start,count);
        Vn=-(k0)^2*((RI).^2-n_imm^2);
    end
    Us=(ifft2((fft2(S.*Vn)).*(AG))*dz); % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
    U=(U);
    %S(:,:,i)=U;
end
evol=U;

    case 'Vol'
        S=zeros(size(V));
for i=1:size(V,3)
    S(:,:,i)=U;
    U=ifft2(propdz.*(fft2(U))); % Prop dz to next layer
    %Us=ifft2(AngGreens(dz).*fftshift(fft2(U.*obj(:,:,i).*dz))); % Scattered field of nth layer
    Us=ifft2((fft2(S(:,:,i).*V(:,:,i))).*(AG))*dz; % Scattered field of nth layer
    %Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
    U=U+Us;
    %S(:,:,i)=U;
end
S(:,:,end+1)=U;
evol=S;

evol=evol(:,:,2:end);

end




% function AG=AngGreensn(Mask,fxx,fyy,lambda,n_imm,prop_phs,z)
% 
% AG= ((-1i.*exp(prop_phs.*z)./(4.*pi.*sqrt((n_imm/lambda)^2-(fxx.^2+fyy.^2))))).*Mask;
% AG(isnan(AG)==1)=0;
% end
end
function Evol=MultiSlice(U_in,RI,lambda,ps,fxx,fyy,n_imm,opt)
% This calculates multiple scattering using the multi-slice method

SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field

fU_current   = fft2(U_in);
prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2));
prop_kernel  = exp(prop_phs * ps);

switch opt
    case 'out'
for layerIdx = 1:RI(2)
    if RI(1) == 0
        RIn = n_imm*ones(size(fxx));
    else
        start=[1 1 layerIdx]; % indicates which layer to read from the data file
        count=[size(fxx) 1]; % Chunk size
        RIn = h5read('sphere.h5','/sphere',start,count);
    end
    Evol  = ifft2(fU_current.* prop_kernel);
    fU_current = fft2(Evol.* exp(1i*2*pi*(RIn)*ps/lambda));
end
        case 'Vol'
Evol=zeros(size(RI));
for layerIdx = 1:size(RI,3)
    Evol(:, :, layerIdx)  = ifft2(fU_current.* prop_kernel);
    fU_current = fft2(Evol(:,:,layerIdx).* exp(1i*2*pi*(RI(:,:,layerIdx))*ps/lambda));
end
end

end
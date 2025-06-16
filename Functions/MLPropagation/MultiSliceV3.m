function Evol=MultiSliceV3(shape,U_in,prop_kernel,inci_flag,k0,dz,n_imm)
% This calculates multiple scattering using the multi-slice method

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [sphere_N(1),  sphere_N(2)];
Nz = sphere_N(3);

fU_current   = fft2(U_in);
%prop_kernel  = P;

for layerIdx = 1:Nz
    if inci_flag == 1
        RIn = n_imm*ones(chunk_size);
    else
        start=[1 1 layerIdx]; % indicates which layer to read from the data file
        count=[chunk_size 1]; % Chunk size
        RIn = h5read('../Medium/sphere.h5','/../Medium/sphere',start,count);
    end
    Evol  = ifft2(fU_current.* prop_kernel);
    fU_current = fft2(Evol.* exp(1i*k0*(RIn)*dz));
end

end
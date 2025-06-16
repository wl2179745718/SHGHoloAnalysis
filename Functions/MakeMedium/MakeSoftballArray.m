function [] = MakeSoftballArray(rj, sigma, n,n_imm, L, X, Y, z,N,delta, k0)

load("../Medium/softballs_parameters.mat")
if (all(softballs_rj==rj,'all')) && (all(softballs_sigma==sigma)) && (all(softballs_n==n)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
    return
end
delete ../Medium/softballs.h5

%[x,y,z] = L2xyz(L,delta);
%[X, Y]=meshgrid(x,y);

Nx = N(1); Ny = N(2); Nz = N(3);

V0 = -(k0)^2.*(n.^2-n_imm^2);

filename='../Medium/softballs.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/softballs', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz
Vn = zeros(size(X));
    for j = 1:size(rj,2)
        Vn= Vn + V0(j) * exp( -((X-rj(1,j)).^2+(Y-rj(2,j)).^2+(z(i)-rj(3,j))^2)/2/sigma(j)^2 );
    end

start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
h5write(filename, '/../Medium/softballs',Vn,start,count);

end

softballs_rj=rj; softballs_sigma=sigma; softballs_n=n; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/softballs_parameters.mat","softballs_rj","softballs_sigma","softballs_n", "Box_L", "Box_delta", "Box_N")

end
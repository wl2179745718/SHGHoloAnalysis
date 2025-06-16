function [] = MakeSoftballHDF5(sigma, n, L, X, Y, z,N,delta, k0)

load("../Medium/softball_parameters.mat")
if (softball_sigma==sigma) && (all(softball_n==n)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
    return
end
delete ../Medium/softball.h5

n_ctr = n(1);
n_imm = n(2);

%[x,y,z] = L2xyz(L,delta);
%[X, Y]=meshgrid(x,y);

Nx = N(1); Ny = N(2); Nz = N(3);

V0 = -(k0)^2*(n_ctr^2-n_imm^2);

filename='../Medium/softball.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/softball', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz

Vn= V0 * exp( -(X.^2+Y.^2+z(i)^2)/2/sigma^2 );
start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
h5write(filename, '/../Medium/softball',Vn,start,count);

end

softball_sigma=sigma; softball_n=n; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/softball_parameters.mat","softball_sigma","softball_n", "Box_L", "Box_delta", "Box_N")

end
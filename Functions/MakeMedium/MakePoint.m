function [] = MakePoint(n, L, x0_index, y0_index, z0_index,N,delta, k0)

load("../Medium/point_parameters.mat")
if (all(sphere_n==n)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
    return
end
delete ../Medium/point.h5

n_sph = n(1);
n_imm = n(2);

%[x,y,z] = L2xyz(L,delta);
%[X, Y]=meshgrid(x,y);

Nx = N(1); Ny = N(2); Nz = N(3);

RI0=n_imm*ones(Nx,Ny);

filename='../Medium/point.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/point', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz

RI=RI0;

if i==z0_index
    RI(x0_index, y0_index)=n_sph;
end

start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
Vn=-(k0)^2*((RI).^2-n_imm^2);
h5write(filename, '/../Medium/point',Vn,start,count);

end

sphere_n=n; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/point_parameters.mat","sphere_n", "Box_L", "Box_delta", "Box_N")

end
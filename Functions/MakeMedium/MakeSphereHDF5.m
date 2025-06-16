function [] = MakeSphereHDF5(rad, n, L, X, Y, z,N,delta, k0)

load("../Medium/sphere_parameters.mat")
if (sphere_rad==rad) && (all(sphere_n==n)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
    return
end
delete ../Medium/sphere.h5

n_sph = n(1);
n_imm = n(2);

%[x,y,z] = L2xyz(L,delta);
%[X, Y]=meshgrid(x,y);

Nx = N(1); Ny = N(2); Nz = N(3);

RI0=n_imm*ones(Nx,Ny);

filename='../Medium/sphere.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/sphere', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz

dielectricSphere=X.^2+Y.^2+z(i)^2<=rad^2;
RI=RI0 + double(dielectricSphere)*(n_sph-n_imm);
Vn=-(k0)^2*((RI).^2-n_imm^2);
start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
h5write(filename, '/../Medium/sphere',Vn,start,count);

end

sphere_rad=rad; sphere_n=n; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/sphere_parameters.mat","sphere_rad","sphere_n", "Box_L", "Box_delta", "Box_N")

end
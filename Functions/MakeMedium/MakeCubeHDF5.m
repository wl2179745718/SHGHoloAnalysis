function [] = MakeCubeHDF5(W, n, L, X, Y, z,N,delta, k0)

%load("../Medium/cube_parameters.mat")
%if (all(cube_W==W)) && (all(cube_n==n)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
%    return
%end
delete ../Medium/cube.h5

n_sph = n(1);
n_imm = n(2);

%[x,y,z] = L2xyz(L,delta);
%[X, Y]=meshgrid(x,y);

Nx = N(1); Ny = N(2); Nz = N(3);

RI0=n_imm*ones(Nx,Ny);

filename='../Medium/cube.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/cube', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz

dielectricSphere=abs(X)<=W(1) & abs(Y)<=W(2) & abs(z(i))<=W(3);
RI=RI0 + double(dielectricSphere)*(n_sph-n_imm);
Vn=-(k0)^2*((RI).^2-n_imm^2);
start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
h5write(filename, '/../Medium/cube',Vn,start,count);

end

cube_W=W; cube_n=n; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/cube_parameters.mat","cube_W","cube_n", "Box_L", "Box_delta", "Box_N")

end
function [] = Make2Sphere(a1, n1, a2, n2, n_imm, theta, R, L, X, Y, z,N,delta, k0)

load("../Medium/TwoSphere_parameters.mat")
if (sphere_a1==a1) && (all(sphere_n1==n1)) && (sphere_a2==a2) && (all(sphere_n2==n2)) && (all(sphere_n_imm==n_imm)) && (sphere_theta==theta) && (all(sphere_R==R)) && (all(Box_L==L)) && (all(Box_delta==delta)) && (all(Box_N==N))
    return
end
delete ../Medium/TwoSphere.h5

Nx = N(1); Ny = N(2); Nz = N(3);

RI0=n_imm*ones(Nx,Ny);

filename='../Medium/TwoSphere.h5';
datasize = [Nx Ny Nz];
h5create(filename, '/../Medium/TwoSphere', datasize, ...
         'Datatype', 'double', 'ChunkSize', [Nx Ny 1] );

for i = 1:Nz

dielectricSphere1=X.^2+Y.^2+z(i)^2<=a1^2;
dielectricSphere2=(X-R*sin(theta)).^2+Y.^2+(z(i)-R*cos(theta))^2<=a2^2;
RI=RI0 + double(dielectricSphere1)*(n1-n_imm)+double(dielectricSphere2)*(n2-n_imm);
Vn=-(k0)^2*((RI).^2-n_imm^2);
start=[1 1 i]; % indicates which layer to read from the data file
count=[Nx Ny 1]; % Chunk size
h5write(filename, '/../Medium/TwoSphere',Vn,start,count);

end

sphere_a1=a1; sphere_n1=n1; sphere_a2=a2; sphere_n2=n2; sphere_n_imm=n_imm; sphere_theta=theta; sphere_R=R; Box_L=L; Box_delta=delta;  Box_N=N;
save("../Medium/TwoSphere_parameters.mat","sphere_a1","sphere_n1","sphere_a2","sphere_n2","sphere_n_imm","sphere_theta","sphere_R", "Box_L", "Box_delta", "Box_N")

end
function [Vi] = sphere_h5read(i,N,n_imm,k0)
% read a layer from the .h5 file

start=[1 1 i]; % indicates which layer to read from the data file
count=[N 1]; % Chunk size
RI = h5read('sphere.h5','/sphere',start,count);
Vi=-(k0)^2*((RI).^2-n_imm^2);

end
function [] = Incident_HDF5(U_inp,N,Pdz)
% Calculate the incident field in the whole space and restore it in a .h5
% file

U=U_inp;

filename='incident.h5';
datasize = [N(1) N(2) N(3)];
h5create(filename, '/real', datasize, ...
         'Datatype', 'double', 'ChunkSize', [N(1) N(2) 1] );
h5create(filename, '/imag', datasize, ...
         'Datatype', 'double', 'ChunkSize', [N(1) N(2) 1] );

for i = 1:N(3)

start=[1 1 i]; % indicates which layer to read from the data file
count=[N(1) N(2) 1]; % Chunk size
h5write(filename, '/real',real(U),start,count);
h5write(filename, '/imag',imag(U),start,count);

U = ifft2(Pdz.*(fft2(U)));

end

end
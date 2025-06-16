function [Uin_i] = Incident_h5read(i,N)
% read a layer from the .h5 file

start=[1 1 i]; % indicates which layer to read from the data file
count=[N 1]; % Chunk size
Uin_i_real = h5read('incident.h5','/real',start,count);
Uin_i_imag = h5read('incident.h5','/imag',start,count);

Uin_i = Uin_i_real + 1i * Uin_i_imag;

end
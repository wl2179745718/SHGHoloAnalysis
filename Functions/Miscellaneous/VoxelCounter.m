function [V_total] = VoxelCounter(shape)
% Count the none-zero voxels of the medium, multiply by Vn

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);
V_total=0;
for i=1:Box_N(3)
    start=[1 1 i]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    V_total=V_total+sum(Vn,'all')*Box_delta(1)*Box_delta(2)*Box_delta(3); % Scattered field of nth layer
    
end

end
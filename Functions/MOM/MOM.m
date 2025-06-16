function U_MOM = MOM(shape, X, Y, box, z, fxin, thetain, lambda, k0)
% Method of moment

box_xy_index = find(abs(X)<box(1)/2 & abs(Y)<box(2)/2);
box_X = X(box_xy_index);
box_Y = Y(box_xy_index);
box_z_index  = find(abs(z)<box(3)/2);
box_z = z(box_z_index);

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);

U_MOM_1d = zeros(1,Box_N(1)*Box_N(2));

for ii = 1:size(box_z_index,2)
    ii
    % Read Vn
    start=[1 1 box_z_index(ii)]; % indicates which layer to read from the data file
    count=[chunk_size 1]; % Chunk size
    Vn = h5read(filename1,filename2,start,count);
    Vn = Vn(box_xy_index);
    % Calculate incident field
    Uin_n = exp(1i*2*pi*fxin*box_X)*exp(1i*2*pi*(box_z(ii)-z(1)+Box_delta(3))*cos(thetain)/lambda);
    % Calculate the Green's function
    parfor jj = 1:Box_N(1)*Box_N(2)
        Rn = sqrt( (X(jj)-box_X).^2 + (Y(jj)-box_Y).^2 + (z(end)-box_z(ii))^2  );
        Gn = -1/4/pi*exp(1i*k0*Rn)./Rn;
        U_MOM_1d(jj) = U_MOM_1d(jj) + sum( Vn.*Uin_n.*Gn );
    end
end

U_MOM = Box_delta(1)*Box_delta(2)*Box_delta(3)*reshape(U_MOM_1d, Box_N(1),  Box_N(2));

end
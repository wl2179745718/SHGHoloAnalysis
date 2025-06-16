function [Us, Us_xz, Us_yz] = MLB4order(shape,U_in,Gdz,G2dz,G3dz,propdz, x0_index, y0_index)

filename = append('../Medium/',shape,'_parameters.mat');
load(filename);
chunk_size = [Box_N(1),  Box_N(2)];

filename1= append('../Medium/',shape,'.h5');
filename2= append('/../Medium/',shape);

% for i=1:4:Box_N(3)
%     %S=U;
%     start=[1 1 i]; % indicates which layer to read from the data file
%     count=[chunk_size 1]; % Chunk size
%     Vn = h5read(filename1,filename2,start,count);
%     X1 = Box_delta(3) * HOlivier(Vn,X0,U,G4dz);
%     U=myifft(propdz.*(myfft(U)));
%     h1 = HOlivier(Vn,X1,U,G3dz);
%     start=[1 1 i+1]; % indicates which layer to read from the data file
%     count=[chunk_size 1]; % Chunk size
%     Vn = h5read(filename1,filename2,start,count);
%     X2 = X1 + Box_delta(3) * h1;
%     U=myifft(propdz.*(myfft(U)));
%     h2 = HOlivier(Vn,X2,U,G2dz);
% 
%     start=[1 1 i+2]; % indicates which layer to read from the data file
%     count=[chunk_size 1]; % Chunk size
%     Vn = h5read(filename1,filename2,start,count);
%     X3 = X2 + Box_delta(3) * h2;
%     U=myifft(propdz.*(myfft(U)));
%     h3 = HOlivier(Vn,X3,U,Gdz);
% 
%     U = 4*Box_delta(3)/3*(2*h1 - h2 + 2*h3) + myifft(propdz.*(myfft(U)));
% 
% end
Ui = U_in;
X0 = zeros(size(U_in));
Us = zeros(size(U_in));
count=[chunk_size 1];

j=1;
for i=1:4:Box_N(3)
    U1 = Us;

    start=[1 1 i];
    Vn = h5read(filename1,filename2,start,count);
    X1 = HOlivier(Vn,X0,Ui+U1,Gdz)*Box_delta(1)*Box_delta(2)*Box_delta(3);
    U1 = myifft(propdz.*(myfft(U1)));
    Ui = myifft(propdz.*(myfft(Ui)));
    start=[1 1 i+1];
    Vn = h5read(filename1,filename2,start,count);
    h1 = HOlivier(Vn,X1,Ui+U1,G3dz);

    X2 = HOlivier(Vn,X1,Ui+U1,Gdz)*Box_delta(1)*Box_delta(2)*Box_delta(3);
    U1 = myifft(propdz.*(myfft(U1)));
    Ui = myifft(propdz.*(myfft(Ui)));
    X1 = myifft(propdz.*(myfft(X1)));
    start=[1 1 i+2];
    Vn = h5read(filename1,filename2,start,count);
    h2 = HOlivier(Vn,X1+X2,Ui+U1,G2dz);

    X3 = HOlivier(Vn,X1+X2,Ui+U1,Gdz)*Box_delta(1)*Box_delta(2)*Box_delta(3);
    U1 = myifft(propdz.*(myfft(U1)));
    Ui = myifft(propdz.*(myfft(Ui)));
    X1 = myifft(propdz.*(myfft(X1)));
    X2 = myifft(propdz.*(myfft(X2)));
    start=[1 1 i+3];
    Vn = h5read(filename1,filename2,start,count);
    h3 = HOlivier(Vn,X1+X2+X3,Ui+U1,Gdz);

    U2 = 4*Box_delta(1)*Box_delta(2)*Box_delta(3)/3*(2*h1 - h2 + 2*h3);
    U1 = myifft(propdz.*(myfft(U1)));
    Ui = myifft(propdz.*(myfft(Ui)));
    Us = U1 + U2;
    Us_yz(j,:) = Us(:,x0_index).';
    Us_xz(j,:) = Us(y0_index,:);
    j=j+1;

end


% for i=1:4:Box_N(3)
%     %S=U;
%     start=[1 1 i]; % indicates which layer to read from the data file
%     Vn = h5read(filename1,filename2,start,count);
%     U1 = Us;
%     X1 = Box_delta(3) * HOlivier(Vn,X0,Ui+U1,G4dz);
%     %U2=myifft(propdz.*(myfft(U2)));
%     Ui=myifft(propdz.*(myfft(Ui)));
%     U1=myifft(propdz.*(myfft(U1)));
% 
%     start=[1 1 i+1]; % indicates which layer to read from the data file
%     Vn = h5read(filename1,filename2,start,count);
%     h1 = HOlivier(Vn,X1,Ui+U1,G3dz);
%     X2 = X1 + Box_delta(3) * h1;
%     Ui=myifft(propdz.*(myfft(Ui)));
%     U1=myifft(propdz.*(myfft(U1)));
% 
%     start=[1 1 i+2]; % indicates which layer to read from the data file
%     Vn = h5read(filename1,filename2,start,count);
%     %U2=myifft(propdz.*(myfft(U2)));
%     h2 = HOlivier(Vn,X2,Ui+U1,G2dz);
%     X3 = X2 + Box_delta(3) * h2;
%     Ui=myifft(propdz.*(myfft(Ui)));
%     U1=myifft(propdz.*(myfft(U1)));
% 
%     start=[1 1 i+3]; % indicates which layer to read from the data file
%     Vn = h5read(filename1,filename2,start,count);
%     %U2=myifft(propdz.*(myfft(U2)));
%     h3 = HOlivier(Vn,X3,Ui+U1,Gdz);
%     Ui=myifft(propdz.*(myfft(Ui)));
%     U1=myifft(propdz.*(myfft(U1)));
% 
%     U2 = 4*Box_delta(3)/3*(2*h1 - h2 + 2*h3);
%     Us = U1 + U2;
% 
% end

end
% Unfinished work: the bug is fixed anyway...

clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 1*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [40*um, 40*um, 20*um];
dr  = [ 80*nm,  80*nm,80*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

Pdz  = Propagator(lambda,FX,FY,dr(3));
propdz = Pdz;
Gdz  = GOlivier(X, Y, dr(1), dr(2),   dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,  dz,dGk,Eps);
AG = Gdz;
G2dz = GOlivier(X, Y, dr(1), dr(2), 2*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,2*dz,dGk,Eps);
G3dz = GOlivier(X, Y, dr(1), dr(2), 3*dr(3), k0);%G_kx_ky(fxx,fyy,n_imm,lambda,3*dz,dGk,Eps);

shape = 'sphere';
nsphere=1.002;%1.2;
rad = 2*um;
figure_window = [-3*rad/um 3*rad/um];

dielectricSphere=X.^2+Y.^2<=rad^2;
RI=n_imm + double(dielectricSphere)*(nsphere-n_imm);
Vn=-(k0)^2*((RI).^2-n_imm^2);

figure
imagesc(x,z,Vn)
colorbar

dz1 = 100*nm;

% H = ifft2((fft2((Us + X).*Vn)).*(AG));

Ui = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);
U1 = zeros(size(Ui));
X0 = zeros(size(Ui));

X1 = myifft((myfft((Ui+U1 + X0).*Vn)).*(AG))*dr(1)*dr(2)*dr(3);

U = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);

S=U;
U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
Us1=myifft((myfft(S.*Vn)).*(AG))*dr(1)*dr(2)*dr(3); % Scattered field of nth layer

norm(X1-Us1,2)/norm(Us1,2)

U1 = myifft(propdz.*(myfft(U1)));
Ui = myifft(propdz.*(myfft(Ui)));
h1 = HOlivier(Vn,X1,Ui+U1,G3dz);

X2 = HOlivier(Vn,X1,Ui+U1,Gdz)*dr(1)*dr(2)*dr(3);
U1 = myifft(propdz.*(myfft(U1)));
Ui = myifft(propdz.*(myfft(Ui)));
X1 = myifft(propdz.*(myfft(X1)));
h2 = HOlivier(Vn,X1+X2,Ui+U1,G2dz);

X3 = HOlivier(Vn,X1+X2,Ui+U1,Gdz)*dr(1)*dr(2)*dr(3);
U1 = myifft(propdz.*(myfft(U1)));
Ui = myifft(propdz.*(myfft(Ui)));
X1 = myifft(propdz.*(myfft(X1)));
X2 = myifft(propdz.*(myfft(X2)));
h3 = HOlivier(Vn,X1+X2+X3,Ui+U1,Gdz);

U2 = 4*dr(1)*dr(2)*dr(3)/3*(2*h1 - h2 + 2*h3);
U1 = myifft(propdz.*(myfft(U1)));
Ui = myifft(propdz.*(myfft(Ui)));
Us = U1 + U2;

figure
subplot(2,2,1)
imagesc(x/um,y/um,abs(X1))
axis square
colorbar
title('|X1|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,2)
imagesc(x/um,y/um,abs(X2))
axis square
colorbar
title('|X2|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,3)
imagesc(x/um,y/um,abs(X3))
axis square
colorbar
title('|X3|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,4)
imagesc(x/um,y/um,abs(U2))
axis square
colorbar
title('|U2|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

U = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);

S=U;
U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
Us1=myifft((myfft(S.*Vn)).*(AG))*dr(1)*dr(2)*dr(3); % Scattered field of nth layer
%Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
U=U+Us1;

S=U;
U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
Us2=myifft((myfft(S.*Vn)).*(AG))*dr(1)*dr(2)*dr(3); % Scattered field of nth layer
%Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
U=U+Us2;

S=U;
U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
Us3=myifft((myfft(S.*Vn)).*(AG))*dr(1)*dr(2)*dr(3); % Scattered field of nth layer
%Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
U=U+Us3;

S=U;
U=myifft(propdz.*(myfft(U))); % Prop dz to next layer
Us4=myifft((myfft(S.*Vn)).*(AG))*dr(1)*dr(2)*dr(3); % Scattered field of nth layer
%Us=conv_fft2(Greens(dz),U.*obj(:,:,i).*dz,'same');
U=U+Us4;

figure
subplot(2,2,1)
imagesc(x/um,y/um,abs(Us1))
axis square
colorbar
title('|Us1|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,2)
imagesc(x/um,y/um,abs(Us2))
axis square
colorbar
title('|Us2|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,3)
imagesc(x/um,y/um,abs(Us3))
axis square
colorbar
title('|Us3|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

subplot(2,2,4)
imagesc(x/um,y/um,abs(Us4))
axis square
colorbar
title('|Us4|')
xlim(figure_window)
ylim(figure_window)
xlabel('x(um)')
ylabel('y(um)')

norm(X1-Us1,2)/norm(Us1,2)
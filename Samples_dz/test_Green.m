clear variables; close all; clc; addpath(genpath('../Functions'));Units

% Box
lambda0 = 1*um; % center wavelength
n_imm = 1;
lambda  = lambda0/n_imm; k0 = 2*pi/lambda;

Box = [10*um, 10*um, 50*um];
dr  = [ 10*nm,  10*nm,10*nm];
[Box,N] = Box_Regularization(Box,dr);

[x, y, z, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);



dz1_array = 10*nm:10*nm:100*nm;

for i = 1:size(dz1_array,2)

dz1=dz1_array(i);
dz2=10*nm;

Pdz  = Propagator(lambda,FX,FY,dz2);

G1 = -exp(1i*k0*sqrt((dz1+dz2)^2+ X.^2 + Y.^2))./sqrt((dz1+dz2)^2+ X.^2 + Y.^2);

G2 = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);

G3 = myifft(Pdz.*(myfft(G2)));

errorG(i) = norm(G1 - G3,2)/norm(G1,2);

end

figure
plot(dz1_array/nm, errorG)
xlabel('dz1(nm)')
ylabel('l2 error of the field')

set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',36)

norm(G1,2)

norm(G3,2)

figure
subplot(1,2,1)
imagesc(abs(G1))
axis square
title('Greens function in r space')
subplot(1,2,2)
FG1 = myfft(G1);
imagesc(fx, fy,abs(FG1).^2)
axis square
title('Greens function in k space')

set(gcf, 'Position', get(0, 'Screensize'));
set(findall(gcf,'-property','FontSize'),'FontSize',36)
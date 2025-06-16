function [] = Diagnose_Greens_funtion_integral_NA_2D_scan_mesh(k0, lambda, NA, meshsize, dz1_array)

Rela_error = zeros(size(dz1_array, 2),size(meshsize,2));
for ii=1:size(dz1_array, 2)
    parfor jj=1:size(meshsize,2)

Box = [320, 320, 5];
dr  = [meshsize(jj), meshsize(jj), 0.01];
dz1 = dz1_array(ii);
figure_window = 20;

[Box,N] = Box_Regularization(Box,dr);

[x, y, z, x0_index, y0_index, z0_index, fx, fy, dfx, dfy, X, Y, FX, FY] = coordinates(N, dr);

Re_fz = real(sqrt(1-FX.^2-FY.^2));
Prop_part = ones(size(Re_fz));
Prop_part(Re_fz==0)=0;

Pdz  = Propagator(lambda,FX,FY,dr(3));

G = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);
%AG= ((-1i.*Propagator(lambda,fxx,fyy,dz)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps))))./dr(1)./dr(2);
GF = myfft(G);

regularization = 2*pi*2/4/pi/dr(1)/dr(2);%abs( GF(x0_index, y0_index) );%(2*pi)^2;% 1/4/pi/dr(1)/dr(2)

GF = GF./regularization;

filter_out = abs( sqrt(FX.^2+FY.^2) ) > NA;
G_NA = GF;
G_NA(filter_out)=0;
Inte_Weyl_NA = dfx * dfy * sum(sum( Re_fz.*abs(G_NA).^2 ));
Inte_Weyl_NA_trapz = trapz(fy,trapz(fx,Re_fz.*abs(G_NA).^2,2));
Analytical = real( 2*pi*(1- sqrt(1-NA^2) ) );
Rela_error(ii,jj) = abs(Inte_Weyl_NA-Analytical)./abs(Analytical).*100;

    end
end

figure
imagesc(meshsize,log10(dz1_array),Rela_error)
set(gca,'YDir','normal')
xlabel('mesh size')
ylabel('log_{10}(distance)')
caxis([0 0.5]);
set(gcf, 'Position', get(0, 'Screensize'));
colorbar
title('% relative error')
function [] = Diagnose_Greens_funtion_integral_NA(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1)

G = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2);
%AG= ((-1i.*Propagator(lambda,fxx,fyy,dz)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps))))./dr(1)./dr(2);
GF = myfft(G);

regularization = 2*pi*2/4/pi/dr(1)/dr(2);%abs( GF(x0_index, y0_index) );%(2*pi)^2;% 1/4/pi/dr(1)/dr(2)

GF = GF./regularization;

index_0 = find(fx>0,1);
index_2 = find(fx>1.2,1);
NA = fx(index_0:index_2);

parfor ii = 1:size(NA,2)
    fprintf(['Diagnosing propagator error accumulation with clipped Greens function: ', num2str([ii]),'/',num2str([size(NA,2)]),'\n'])
    filter_out = abs( sqrt(FX.^2+FY.^2) ) > NA(ii);
    G_NA = GF;
    G_NA(filter_out)=0;
    Inte_Weyl_NA(ii) = dfx * dfy * sum(sum( Re_fz.*abs(G_NA).^2 ));
    Inte_Weyl_NA_trapz(ii) = trapz(fy,trapz(fx,Re_fz.*abs(G_NA).^2,2));
    Analytical(ii) = real( 2*pi*(1- sqrt(1-NA(ii)^2) ) );
end

figure
subplot(1,2,1)
imagesc(fx,fy,abs(GF))

subplot(1,2,2)
imagesc(X(1,:),Y(:,1),abs(G))


figure
yyaxis left
hold on
plot(NA, abs(Inte_Weyl_NA))
plot(NA, abs(Analytical))
plot(NA, abs(Inte_Weyl_NA-Analytical))
ylim([0 8]);
xlabel('NA')
yyaxis right
plot(NA, abs(Inte_Weyl_NA-Analytical)./abs(Analytical).*100)
ylabel('% error')
ylim([0 10]);
%ylabel('\int \frac{1}{\gamma}|G|^2 d^2 k_\perp')
legend('$\int Re(\gamma)|G|^2 d^2 k_\perp$','$Re( 2\pi(1-\sqrt{1-NA^2}) )$','difference','Interpreter','latex');
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Integration by Summation')

figure
yyaxis left
hold on
plot(NA, abs(Inte_Weyl_NA_trapz))
plot(NA, abs(Analytical))
plot(NA, abs(Inte_Weyl_NA_trapz-Analytical))
ylim([0 8]);
xlabel('NA')
yyaxis right
plot(NA, abs(Inte_Weyl_NA-Analytical)./abs(Analytical).*100)
ylabel('% error')
ylim([0 10]);
%ylabel('\int \frac{1}{\gamma}|G|^2 d^2 k_\perp')
legend('$\int Re(\gamma)|G|^2 d^2 k_\perp$','$Re( 2\pi(1-\sqrt{1-NA^2}) )$','difference','Interpreter','latex');
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Integration by 2D trapz')


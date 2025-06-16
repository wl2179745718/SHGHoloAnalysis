function [] = Diagnose_Greens_funtion_integral(X, Y, z_positive, x0_index, y0_index, FX, FY, fx, fy,Re_fz, dfx, dfy, Prop_part, dr, Pdz, k0, lambda, dz1)

G2 = -exp(1i*k0*sqrt((dz1)^2+ X.^2 + Y.^2))./sqrt((dz1)^2+ X.^2 + Y.^2); 
G2F = myfft(G2);
barMax = max(max( log(abs(G2F)) ));
barMin = log(abs(G2F(round(end/2),round(end/2))))-3;
%G3_0 = G2F;

regularization = abs( G2F(x0_index, y0_index) );

z_positive = z_positive - dr(3) + dz1;

for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation with clipped Greens function: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2); 
G1F = myfft(G1);

Inte_G1F(ii) = dfx * dfy / regularization * sum(sum( Prop_part.*abs(G1F).^2 ));
Inte_Weyl_G1F(ii) = dfx * dfy / regularization * sum(sum( Re_fz.*abs(G1F).^2 ));


errorG(ii) = norm(G1 - G2,2)/norm(G1,2);

G2F = Pdz.*(myfft(G2));
%G3k = Propagator(lambda,FX,FY,z_positive(ii)-dz1+dr(3)).*(G3_0);

Inte_G2F(ii) = dfx * dfy / regularization * sum(sum( Prop_part.*abs(G2F).^2 ));
Inte_Weyl_G2F(ii) = dfx * dfy / regularization * sum(sum( Re_fz.*abs(G2F).^2 ));

G2 = myifft(G2F); 
%G3 = myifft(G3k);

%errorG2G3(ii) = norm(G2-G3,2)/norm(G2,2);

%G2k_clip = G2k;
%G2k_clip( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;
%G3k_clip = G3k;
%G3k_clip( lambda^2* (FX.^2+FY.^2) >=0.8 )=0;
%G2_clip = myifft(G2k_clip); 
%G3_clip = myifft(G3k_clip);

%errorG2G3_clip(ii) = norm(G2_clip-G3_clip,2)/norm(G2_clip,2);

end

figure('Name','Diagnose: Integration of G1')
subplot(2,1,1)
plot(z_positive/lambda,Inte_G1F)
text(z_positive(end)/lambda/2, 0.2*max(errorG)+0.8*min(errorG), ['dz1 = ', num2str(dz1/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('\int G1 df_x df_y')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose Greens function integration.png");
imshow(A,map)

figure('Name','Diagnose: Weyl Integration of G1')
subplot(2,1,1)
plot(z_positive/lambda,Inte_Weyl_G1F)
text(z_positive(end)/lambda/2, 0.2*max(errorG)+0.8*min(errorG), ['dz1 = ', num2str(dz1/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('\int G1 df_x df_y')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose Greens function integration.png");
imshow(A,map)

figure('Name','Diagnose: Integration of G2')
subplot(2,1,1)
plot(z_positive/lambda,Inte_G2F)
text(z_positive(end)/lambda/2, 0.2*max(errorG)+0.8*min(errorG), ['dz1 = ', num2str(dz1/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('\int G1 df_x df_y')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose Greens function integration.png");
imshow(A,map)

figure('Name','Diagnose: Weyl Integration of G2')
subplot(2,1,1)
plot(z_positive/lambda,Inte_Weyl_G2F)
text(z_positive(end)/lambda/2, 0.2*max(errorG)+0.8*min(errorG), ['dz1 = ', num2str(dz1/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('\int G1 df_x df_y')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Diagnose Greens function integration.png");
imshow(A,map)

end
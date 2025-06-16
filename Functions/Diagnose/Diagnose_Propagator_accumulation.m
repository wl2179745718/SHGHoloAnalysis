function [] = Diagnose_Propagator_accumulation(X, Y, z_positive, dr, Pdz, k0, lambda)

G2 = -exp(1i*k0*sqrt((dr(3))^2+ X.^2 + Y.^2))./sqrt((dr(3))^2+ X.^2 + Y.^2); 
for ii = 1:size(z_positive,2)

fprintf(['Diagnosing propagator error accumulation: ', num2str([ii]),'/',num2str([size(z_positive,2)]),'\n'])

G1 = -exp(1i*k0*sqrt((z_positive(ii))^2+ X.^2 + Y.^2))./sqrt((z_positive(ii))^2+ X.^2 + Y.^2); 

errorG(ii) = norm(G1 - G2,2)/norm(G1,2);

G2 = myifft(Pdz.*(myfft(G2))); 

end

figure('Name','Diagnose: propagator error accumulation')
subplot(2,1,1)
plot(z_positive/lambda,errorG)
text(z_positive(end)/lambda/2, 0.9*max(errorG), ['dz = ', num2str(dr(3)/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Propagation error accumulation.png");
imshow(A,map)

end
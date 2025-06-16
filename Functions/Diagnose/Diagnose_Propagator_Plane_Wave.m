function [] = Diagnose_Propagator_Plane_Wave(uin,z,dr,thetain,lambda,Pdz)

u2 = uin;
for ii=1:size(z,2)
    fprintf(['Diagnosing propagator error with plane wave: ', num2str([ii]),'/',num2str([size(z,2)]),'\n'])
    u1 = uin*exp(1i*2*pi*(z(ii)-z(1)+dr(3))*cos(thetain)/lambda);
    u2 = myifft(Pdz.*(myfft(u2))); 
    %abs(u2(x0_index, y0_index))
    errorU(ii) = norm(u1 - u2,2)/norm(u1,2);
end

figure('Name','Diagnose: propagator with plane wave')
subplot(2,1,1)
plot(z/lambda,errorU)
%text(z(end)/lambda/2, 0.9*max(errorU), ['dz = ', num2str(dr(3)/lambda),' \lambda'], 'FontSize', 30); 
xlabel('z(\lambda)')
ylabel('error')
set(gcf, 'Position', get(0, 'Screensize'));
subplot(2,1,2)
[A,map] = imread("Plane wave propagation.png");
imshow(A,map)

end
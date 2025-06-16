function [] = Diagnose_Green(Green, fx, fy, lambda)
% The Green's funciton is in k-space

src = zeros(size(Green));
src(round(end/2),round(end/2))=1;

figure('Name','Diagnose: error in the Greens function')
subplot(1,3,1)
imagesc(fx.*lambda, fy.*lambda,abs(Green))
axis square
title('|G| in k space')
colorbar

subplot(1,3,2)
plot(fx.*lambda,abs(Green(round(end/2),:)))
xlim([-2 2])
title('|G| vs. kx')

subplot(1,3,3)
plot(fx.*lambda,angle(Green(round(end/2),:)))
xlim([-2 2])
title('\angle G vs. kx')

set(gcf, 'Position', get(0, 'Screensize'));

end
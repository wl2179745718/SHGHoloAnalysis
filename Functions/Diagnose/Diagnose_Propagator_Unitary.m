function [] = Diagnose_Propagator_Unitary(lambda,FX,FY,dr,x0_index,fy)

Pdz  = Propagator(lambda,FX,FY,dr(3));
P2dz  = Propagator(lambda,FX,FY,2*dr(3));
P3dz  = Propagator(lambda,FX,FY,3*dr(3));

Plambda  = Propagator(lambda,FX,FY,lambda);
P2lambda  = Propagator(lambda,FX,FY,2*lambda);
P4lambda  = Propagator(lambda,FX,FY,4*lambda);

Pdz_1d = Pdz(:,x0_index);
P2dz_1d = P2dz(:,x0_index);
P3dz_1d = P3dz(:,x0_index);

Plambda_1d = Plambda(:,x0_index);
P2lambda_1d = P2lambda(:,x0_index);
P4lambda_1d = P4lambda(:,x0_index);

figure('Name','Diagnose: Propagator of different distances')
subplot(2,3,1)
plot(lambda*fy,abs(Pdz_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{\Delta z}|')
xlabel('\lambda f_y')

subplot(2,3,2)
plot(lambda*fy,abs(P2dz_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{2\Delta z}|')
xlabel('\lambda f_y')

subplot(2,3,3)
plot(lambda*fy,abs(P3dz_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{3\Delta z}|')
xlabel('\lambda f_y')

subplot(2,3,4)
plot(lambda*fy,abs(Plambda_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{\lambda}|')
xlabel('\lambda f_y')

subplot(2,3,5)
plot(lambda*fy,abs(P2lambda_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{2\lambda}|')
xlabel('\lambda f_y')

subplot(2,3,6)
plot(lambda*fy,abs(P4lambda_1d))
xlim([-2 2])
ylim([0 1])
title('|P_{4\lambda}|')
xlabel('\lambda f_y')

set(gcf, 'Position', get(0, 'Screensize'));

end
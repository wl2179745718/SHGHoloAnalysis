clear variables; close all; clc; addpath(genpath('../Functions'));
%https://core.ac.uk/download/pdf/30889539.pdf

lambda = 1;
k = 2*pi/lambda;
m = 1:0.01:2;

%a = 0:0.1:40;
x = 3; %k*a;
rad = x./k;
rho = 2*x*(m-1);

n_max = 40;
n = 1:n_max;

%S0_n = (2*n+1).*(an+bn);

%Q_ext_Mie = 4/k^2./a^2*Re(S0);

for ii = 1:size(m,2)
[ an, bn ] = sph_exp_coef_mie( n_max, lambda, m(ii), rad );
S0_n = (2*n+1).*(an+bn); %Eq 3.6
S0(ii) = 1/2* sum(S0_n);
end
Q_ext_Mie = 4/k^2./rad.^2.*real(S0);

Q_ext_apr = 2 - 4./rho.*sin(rho)+4./rho.^2.*(1-cos(rho));

figure(100)
plot(m,Q_ext_apr)
hold on
plot(m,Q_ext_Mie)
title('Replicate Fig 3.3 and compare with simu results')

%x_simu = 0:2:20; %k*a;
%rad_simu = x_simu./k;
m_simu = 1:0.04:2;
k0 = 1/lambda;

ord = 1;
for jj = 1:size(m_simu,2)

ps          = .1;                 % pixel size (x,y,z) in object space (microns)
n_imm       = 1;                % refractive index of immersion media
nsphere= m(jj);%1.2;
n=[nsphere,n_imm];
N           = [2^9, 2^9, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.4;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
%rad = rad_simu(jj);

RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);

%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.25;

U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));


E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,'Vol');
%E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
%E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,'Vol');%.017
%E_MSL_3d=MultiSlice(U_inp,RI,lambda,ps,fxx,fyy,n_imm);

V_inc = 0*V;
RI_inc = n_imm * ones(size(RI));
E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,'Vol');
%E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
%E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,'Vol');%.017
%E_MSL_3d_inc=MultiSlice(U_inp,RI_inc,lambda,ps,fxx,fyy,n_imm);
U_inp_end_MSR = E_MSR_3d_inc(:,:,end);
%U_inp_end_MLR = E_MLR_3d_inc(:,:,end);
%U_inp_end_MLB = E_MLB_3d_inc(:,:,end);
%U_inp_end_MSL = E_MSL_3d_inc(:,:,end);


[E_tot_MSR,E_sca_MSR] = tot2sca(E_MSR_3d,U_inp_end_MSR,'Vol');
%[E_tot_MLR,E_sca_MLR] = tot2sca(E_MLR_3d,U_inp_end_MLR);
%[E_tot_MLB,E_sca_MLB] = tot2sca(E_MLB_3d,U_inp_end_MLB);
%[E_tot_MSL,E_sca_MSL] = tot2sca(E_MSL_3d,U_inp_end_MSL);

S0_simu(jj) = abs(E_sca_MSR(end/2,end/2))^2;

end
%Normalize
Q_ext_simu = 4/k^2./rad.^2.*real(S0_simu);

Q_ext_simu = 50*Q_ext_simu./max(Q_ext_simu);


figure(100)
plot(m_simu,Q_ext_simu,'-o');

figure
imagesc(E_sca_MSR)
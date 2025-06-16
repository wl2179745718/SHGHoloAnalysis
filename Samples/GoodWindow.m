clear variables; close all; clc; addpath(genpath('../Functions'));
% unit = m

N_array = [2^8, 2^9, 3*2^8, 2^10, 3*2^9];
error_MLB_l2 = zeros(size(N_array));
error_MLB_max= zeros(size(N_array));
error_MLR_l2 = zeros(size(N_array));
error_MLR_max= zeros(size(N_array));
error_MSR_l2 = zeros(size(N_array));
error_MSR_max= zeros(size(N_array));
error_MSL_l2 = zeros(size(N_array));
error_MSL_max= zeros(size(N_array));

for ii = 1:size(N_array,2)

ps          = .1;       % pixel size (x,y,z) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
nsphere=1.01;%1.2;
n=[nsphere,n_imm];
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [N_array(ii)-1, N_array(ii)-1, 2^9];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.6;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

opt = 'Vol';

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 1;
rad = 2;

switch opt
    case 'Vol'
RI = MakeSphereInRandMed(rad, n, L, delta);
V=-(k0)^2*((RI).^2-n_imm^2);
V_inc = 0*V;
RI_inc = n_imm * ones(size(RI));
    case 'out'
RI = [1 N(3)];
V  = RI;
MakeSphereHDF5(rad, n, L, delta);
RI_inc = [0 N(3)];
V_inc = RI_inc;
end
%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.05;


%start=[1 1 280]; % indicates which layer to read from the data file
%count=[size(fxx) 1]; % Chunk size
%RI0 = h5read('sphere.h5','/sphere',start,count);

%figure
%imagesc(x,y,RI0)
%axis square
%colorbar


U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

ord = 1;

E_MSR_3d=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,ord,Eps,dGk,opt);
E_MLR_3d=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,opt);%.017
E_MLB_3d=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V,U_inp,Eps,dGk,opt);%.017
E_MSL_3d=MultiSlice(U_inp,RI,lambda,ps,fxx,fyy,n_imm,opt);

E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,opt);
E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,opt);%.017
E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,opt);%.017
E_MSL_3d_inc=MultiSlice(U_inp,RI_inc,lambda,ps,fxx,fyy,n_imm,opt);

switch opt
    case 'Vol'
U_inp_end_MSR = E_MSR_3d_inc(:,:,end);
U_inp_end_MLR = E_MLR_3d_inc(:,:,end);
U_inp_end_MLB = E_MLB_3d_inc(:,:,end);
U_inp_end_MSL = E_MSL_3d_inc(:,:,end);
    case 'out'
U_inp_end_MSR = E_MSR_3d_inc;
U_inp_end_MLR = E_MLR_3d_inc;
U_inp_end_MLB = E_MLB_3d_inc;
U_inp_end_MSL = E_MSL_3d_inc;
end

%E_MSR_xz = squeeze(E_MSR_3d(round(end/2),:,:));
%E_MLR_xz = squeeze(E_MLR_3d(round(end/2),:,:));
%E_MLB_xz = squeeze(E_MLB_3d(round(end/2),:,:));
%E_MSL_xz = squeeze(E_MSL_3d(round(end/2),:,:));


[E_tot_MSR,E_sca_MSR] = tot2sca(E_MSR_3d,U_inp_end_MSR,opt);
[E_tot_MLR,E_sca_MLR] = tot2sca(E_MLR_3d,U_inp_end_MLR,opt);
[E_tot_MLB,E_sca_MLB] = tot2sca(E_MLB_3d,U_inp_end_MLB,opt);
[E_tot_MSL,E_sca_MSL] = tot2sca(E_MSL_3d,U_inp_end_MSL,opt);

E_sca_MLB = E_sca_MLB./max(max(E_sca_MLB));
E_sca_MLR = E_sca_MLR./max(max(E_sca_MLR));
E_sca_MSR = E_sca_MSR./max(max(E_sca_MSR));
E_sca_MSL = E_sca_MSL./max(max(E_sca_MSL));

[x_sca_MLB, y_sca_MLB] = center(X,Y,E_sca_MLB);
[x_sca_MLR, y_sca_MLR] = center(X,Y,E_sca_MLR);
[x_sca_MSR, y_sca_MSR] = center(X,Y,E_sca_MSR);
[x_sca_MSL, y_sca_MSL] = center(X,Y,E_sca_MSL);


[E_sca_Mie] = Mie_plane(X,Y,z(end),phi,k,rad,c/lambda,n_imm^2,1,nsphere^2,1,40);
[x_sca_Mie, y_sca_Mie] = center(X,Y,E_sca_Mie);


figure
subplot(2,2,1)
hold on
plot(unwrap(angle(E_sca_MLB(round(end/2),:))))
plot(unwrap(angle(E_sca_Mie(round(end/2),:))))
subplot(2,2,2)
plot(unwrap(angle(E_sca_MLB(:,round(end/2)))))
plot(unwrap(angle(E_sca_Mie(:,round(end/2)))),'o')
subplot(2,2,3)
plot(abs(E_sca_MLB(:,round(end/2))))
plot(abs(E_sca_Mie(:,round(end/2))),'o')
subplot(2,2,4)
plot(abs(E_sca_MLB(round(end/2),:)))
plot(abs(E_sca_Mie(round(end/2),:)),'o')

[~,m_index]=max(abs(E_sca_MLB(round(end/2),:)));
E_sca_MLB(round(end/2),m_index)
x(m_index)

[m_value,m_index]=max(abs(E_sca_Mie(round(end/2),:)));
E_sca_Mie(round(end/2),m_index)
x(m_index)

%corr_m = corr(E_sca_MLB(round(end/2),:), E_sca_Mie(round(end/2),:));
crr = xcorr2(E_sca_MLB,E_sca_Mie);
[ssr,snd] = max(crr(:));
[i,j] = ind2sub(size(crr),snd)


thrsh = 0.2;

%cond = abs(E_sca_Mie) > thrsh & abs(E_sca_MSL) > thrsh;

%E_Mie_corr = E_sca_Mie.*cond;
%E_MSL_corr = E_sca_MSL.*cond;

%fun = @(x)norm(abs((x(1)+1i*x(2))*E_Mie_corr - E_MSL_corr));
%x0 = [1,0];
%x_factor = fminsearch(fun,x0);

%E_Mie_MSL = (x_factor(1)+1i*x_factor(2))*E_sca_Mie;

E_sca_MLB = normalize(E_sca_MLB, E_sca_Mie, thrsh);
E_sca_MLR = normalize(E_sca_MLR, E_sca_Mie, thrsh);
E_sca_MSR = normalize(E_sca_MSR, E_sca_Mie, thrsh);
E_sca_MSL = normalize(E_sca_MSL, E_sca_Mie, thrsh);

error_MLB_l2(ii) = norm(abs(E_sca_MLB-E_sca_Mie),2)./norm(abs(E_sca_Mie),2);
error_MLB_max(ii)= max(max(abs(E_sca_MLB-E_sca_Mie)));
error_MLR_l2(ii) = norm(abs(E_sca_MLR-E_sca_Mie),2)./norm(abs(E_sca_Mie),2);
error_MLR_max(ii)= max(max(abs(E_sca_MLR-E_sca_Mie)));
error_MSR_l2(ii) = norm(abs(E_sca_MSR-E_sca_Mie),2)./norm(abs(E_sca_Mie),2);
error_MSR_max(ii)= max(max(abs(E_sca_MSR-E_sca_Mie)));
error_MSL_l2(ii) = norm(abs(E_sca_MSL-E_sca_Mie),2)./norm(abs(E_sca_Mie),2);
error_MSL_max(ii)= max(max(abs(E_sca_MSL-E_sca_Mie)));

%delete sphere.h5

end

figure
subplot(2,4,1)
plot(ps*N_array,error_MLB_l2,'o');
title('MLB')
subplot(2,4,2)
plot(ps*N_array,error_MLR_l2,'o');
title('MLR')
subplot(2,4,3)
plot(ps*N_array,error_MSR_l2,'o');
title('MSR')
subplot(2,4,4)
plot(ps*N_array,error_MSL_l2,'o');
title('MSL')
subplot(2,4,5)
plot(ps*N_array,error_MLB_max,'o');
subplot(2,4,6)
plot(ps*N_array,error_MLR_max,'o');
subplot(2,4,7)
plot(ps*N_array,error_MSR_max,'o');
subplot(2,4,8)
plot(ps*N_array,error_MSL_max,'o');
saveas(gcf,'Goodwindow.png')

clear variables; close all; clc; addpath(genpath('../Functions'));
% unit = m
ps          = .1;       % pixel size (x,y,z) in object space (microns)
lambda      =  1;       % central wavelength (microns)  
n_imm       =  1;       % refractive index of immersion media
nsphere=1.02;%1.2;
n=[nsphere,n_imm];
c=299792458;
k0=(2*pi)/lambda;
k=k0*n_imm;
N           = [3, 3, 2];                  % lateral pixel dimension 
L = ps*N;
delta = [ps, ps, ps];
deltaf = 1./L;
NA_in = 0.0;
NA_in = round(NA_in/deltaf(1)*n_imm/lambda)*deltaf(1)/n_imm*lambda;
phi       = asin(NA_in);              

opt = 'Vol';

[x,y,z] = L2xyz(L,delta);
[X,Y]=meshgrid(x,y);
[fx,fy] = L2fxfy(L,delta);
[fxx,fyy]   = meshgrid(fx,fx);      % 2D grid in fx/fy

dGk = 0;
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
V_inc = RI;
end
%z_inc = z(1)-ps;
Eps=n_imm^2/lambda^2*0.25;


U_inp=exp(1i*k*sin(phi)*X);%ones(N(1),N(2));

ord = 1;


%E_MSR_3d_inc=MultiSlabRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,ord,Eps,dGk,opt);
%E_MLR_3d_inc=MultiLayerRytovv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,opt);%.017
%E_MLB_3d_inc=MultiLayerBornv2(fxx,fyy,lambda,n_imm,ps,V_inc,U_inp,Eps,dGk,opt);%.017
%E_MSL_3d_inc=MultiSlice(U_inp,RI_inc,lambda,ps,fxx,fyy,n_imm,opt);

%E_MLB_3d(:,:,end)

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


opt = 'out';

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
U_inp_end_MSR2d = E_MSR_3d_inc;
U_inp_end_MLR2d = E_MLR_3d_inc;
U_inp_end_MLB2d = E_MLB_3d_inc;
U_inp_end_MSL2d = E_MSL_3d_inc;
end

%E_MSR_xz = squeeze(E_MSR_3d(round(end/2),:,:));
%E_MLR_xz = squeeze(E_MLR_3d(round(end/2),:,:));
%E_MLB_xz = squeeze(E_MLB_3d(round(end/2),:,:));
%E_MSL_xz = squeeze(E_MSL_3d(round(end/2),:,:));


norm(U_inp_end_MSR-U_inp_end_MSR2d)
norm(U_inp_end_MLR-U_inp_end_MLR2d)
norm(U_inp_end_MLB-U_inp_end_MLB2d)
norm(U_inp_end_MSL-U_inp_end_MSL2d)
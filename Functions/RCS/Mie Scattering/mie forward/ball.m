%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file to simulate holograms from sparsely distributed particles in 3D
% The scattering from particles are generated according to Mie theory.
% Matlab implementation of Mie scattering was taken from work from 
% http://omlc.org/software/mie/
% by Christian Maetzler in 2002.
%
% reference:
% Lei Tian, Compressive Phase retrieval, Chapter 2 of PhD Thesis, MIT, 2013. 
% Wensheng Chen, et. al, 
% "Empirical concentration bounds for compressive holographic bubble 
% imaging based on a Mie scattering model," Opt. Express 23, 4715-4725 (2015)
%
% The program allows particles with random diameter distributions, with
% different index refractions, and distributed in 3D random locations.
%
% last modified on 10/07/2015
% by Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
path([pwd,filesep,'MieFunctions'],path);
if ~exist("my_parpool", 'var')
    my_parpool = parpool('threads');
end

%% system parameters: NA needs to be large enough
% size of hologram
N = 4000;
% pixel size
dpix = 1; % length in px

d = 24; % 24px=6*lambda


%% parameters of the particles
% medium
n1 = 1.;   % index of refraction of water
% index of the particle
n2 = 1.02;

lambda = 4; % in px
k = 2*pi/lambda*n1;

% delta_x = lambda/NA % lateral resolution
% DOF = lambda/NA^2 % depth of focus
% all_ill=[0.8906666666666667, 0.0; 0.6293333333333333, 0.6293333333333333; 0.0, 0.8906666666666667; -0.6293333333333333, 0.6293333333333333; -0.8906666666666667, 0.0; -0.6293333333333333, -0.6293333333333333; 0.0, -0.8906666666666667; 0.6293333333333333, -0.6293333333333333];
all_ill=[0 0];
all_ill=all_ill(1,:);

z = 500;
LEDs=size(all_ill, 1);

%Esave = zeros(Hsize,Hsize,length(z));
%% Mie solver
% assuming x-polarized illumination
% generating total field0

pupil=load("../4k_0.9pupil.mat").data;
mask=load("../4k_mask.mat").data;
prop=load("../4k_prop-500.mat").data;

% mask=im2double(imread("../3k_crop.tiff"));
% pupil=load("../3k_0.9pupil.mat").data;
% prop=load("../prop_3000_-1000.mat").data;


Holo=zeros(N,N,LEDs);


tic
for i=1:LEDs
    i
    c_ab = all_ill(i, :);
    [E, inv_phase] = Mie_x_radiation2(n1, n2, d, lambda, N, dpix, z, [-0.5,-0.5], [c_ab,sqrt(1-c_ab*c_ab')], 5);
    % [E, inv_phase] = Mie_x_radiation2(n1, n2, d, lambda, N, dpix, z, [-0.5,-0.5], [0,0,1], 6);
    Etot = ifft2(fft2(E(:,:,1).*mask).*pupil).*inv_phase;
    % generating hologram
    Holo(:,:,i) = 1+2*real(Etot)+abs(Etot).^2;
    % record totol field for comparison.
    % Field(:,:,i) = Etot;
end
toc
figure; imagesc(squeeze(Holo(:,:,1))); axis image; caxis([0.9 1.1])
% 
% fn = ['MieBall',num2str(N)];
% save(fn, 'Holo', 'd', 'dpix', 'lambda','illum', 'n1', 'n2');

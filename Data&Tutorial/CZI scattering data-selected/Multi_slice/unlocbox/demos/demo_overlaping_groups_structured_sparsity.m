%DEMO_OVERLAPING_GROUPS_STRUCTURED_SPARSITY Demonstration for the use of overlaping structured sparsity
%
%   Here we solve a sound declipping problem. The problem can be
%   expressed as this 
%
%        argmin_x   || x ||_21 such that A G^* x = b
%
%   where b is the signal at the non clipped part,  A an operator
%   representing the mask selecting the non clipped part of the signal and
%   G^ is the Gabor synthesis operation
%
%   Here the general assumption is that the signal is sparse in the Gabor
%   domain!
%
%   Warning! Note that this demo requires the LTFAT toolbox to work.
%
%   We set 
%
%    f_1(x)=||x||_{21}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||z||_21
%
%     The groups are defined like this
%
%     For a 2 by 8 spectrogram 
%     11112222        21111222        22111122        22211112       
%     33334444        43333444        44333344        44433334
%
%    f_2(x)=||Ax-b||_2^2
%     We define the gradient as: 
%
%        grad_f(x) = 2 * G A^*( A G^*x - b )
%
%   Results
%   -------
%
%   Figure 1: Original spectrogram
%
%      This figure shows the original spectrogram.
%
%   Figure 2: Spectrogram of the depleted sound
%
%      This figure shows the spectrogram after the loss of the sample (We loos 75% of the samples.)
%
%   Figure 3: Spectrogram of the reconstructed sound
%
%      This figure shows the spectrogram of the reconstructed sound thanks to the algorithm.
%
%   References:
%     H. Raguet, J. Fadili, and G. Peyré. Generalized forward-backward
%     splitting. arXiv preprint arXiv:1108.4404, 2011.
%     
%     K. Siedenburg, M. Kowalski, and M. Dorfler. Audio declipping with
%     social sparsity. In Acoustics, Speech and Signal Processing (ICASSP),
%     2014 IEEE International Conference on, pages 1577--1581. IEEE, 2014.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/demos/demo_overlaping_groups_structured_sparsity.php

% Copyright (C) 2012-2013 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.6.2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Author: Nathanael Perraudin
% Date: January 8 2015
%



%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();
ltfatstart(); % start the ltfat toolbox

verbose = 2;    % verbosity level
writefile=0;    % writting wav sound

%% Defining the problem

% Original sound
[sound_original, Fs]=gspi();



%%
length_sig=length(sound_original); % Put a small number here if you want to proceed only a part a of the signal
sound_part=sound_original(1:length_sig);

% In oder to write the depleted sound somewhere
if writefile
    wavwrite(sound_part,Fs,'original.wav');
end

tmax = 0.08;
tmin = -0.3;
Mask = 1-(sound_part>tmax) - (sound_part<tmin);

% Depleted sound
sound_depleted = sound_part;
sound_depleted(sound_part>tmax) = tmax;

sound_depleted(sound_part<tmin) = tmin;
if writefile
    wavwrite(sound_depleted,Fs,'depleted.wav');
end

%% Setting proximal operators

tau = 1e-2; % regularization parameter for the problem

% select a gabor frame for a real signal with a Gaussian window
a=64; % size of the shift in time
M=256;% number of frequencies
F=frametight(frame('dgtreal','gauss',a,M));
% Get the framebounds
GB = M/a;

% Define the Frame operators
Psi = @(x) frana(F,x);
Psit = @(x) frsyn(F,x);

f2.prox = @(x,T) Psi(Psit(x) .* ( 1 -  Mask )+ Mask.* sound_depleted);
f2.eval = @(x) eps;

% setting the function f1 (l1 norm of the Gabor transform)
% set parameters
param_l12.verbose = verbose - 1;



% The groups are made like this
%
%   For a 2 by 8 spectrogram 
%   11112222        21111222        22111122        22211112       
%   33334444        43333444        44333344        44433334

% Create the group (this is a bad code)
    % -------------------------------------------- %
    lg = 4; % length of the group;
    
    xin = Psi(sound_depleted);
    xin_im = framecoef2native(F,xin);
    
    K=size(xin_im,1);
    L=size(xin_im,2);
    
    g_d1=zeros(K*L,1);
    
    indice = 1:length(g_d1);
    
    indice_mat =reshape(indice,L,K)';
    sgd = floor(L/lg)*lg * K;
    g_d = zeros(lg,sgd);
    g_t = lg*ones(lg,sgd/lg);
    
    for ii=1:lg
       for jj=1 : floor(L/lg)
           for ll = 1:K
                g_d( ii, (1:4) + (jj-1)*lg + floor(L/lg)*lg * (ll-1)) = ...
                indice_mat(ll, mod((1:4) + (jj-1)*lg+ii -1, L)+1 );
           end            
       end
    end




    
    
    % -------------------------------------------- %
  
param_l12.g_t = g_t;
param_l12.g_d = g_d;


f1.prox=@(x, T) prox_l21(x, T*tau, param_l12);
f1.eval=@(x) tau*norm_l21(x,g_d, g_t);   



%% solving the problem


% setting different parameters  for the simulation
param.verbose = verbose; % display parameter
param.maxit = 100; % maximum iteration
param.tol = 1e-6; % tolerance to stop iterating
%param.method = 'FISTA'; % desired method for solving the problem

%sol=Psit(forward_backward(xin,f1,f2,param));

param.do_ts = @(x) log_decreasing_ts(x, 10, 0.1, 80);
sol=Psit(solvep(Psi(sound_part),{f1,f2},param));




%% Evaluate the result
snr_in = snr(sound_part,sound_depleted);
snr_fin = snr(sound_part,sol);


fprintf('The SNR of the initial signal is %g dB \n',snr_in);
fprintf('The SNR of the recovered (FB) signal is %g dB \n',snr_fin);



% In order to write the restored sound somewhere
if writefile
    wavwrite(sol,Fs,'restored.wav');
end
%%
dr=90;

figure(1);
plotframe(F,Psi(sound_part),Fs,dr);
title('Gabor transform of the original sound');

figure(2);
plotframe(F,Psi(sound_depleted),Fs,dr);
title('Gabor transform of the depleted sound');

figure(3);
plotframe(F,Psi(sol),Fs,dr);
title('Gabor transform of the reconstructed sound');



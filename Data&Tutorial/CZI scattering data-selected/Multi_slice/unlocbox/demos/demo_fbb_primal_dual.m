%DEMO_FBB_PRIMAL_DUAL  Example of use of the forward backward based primal dual  solver
%
%   We present an example of the the forward backward based primal dual
%   solver through an image de-noising, in-painting problem. We express the
%   problem in the following way
%
%        argmin ||A(x-b)||^2 + lambda*||x||_TV + tau*||Wx||_1
%
%   Where b is the degraded image, W the wavelet transform and A a
%   linear operator performing the masking operation. This operator set to
%   0 all unknown pixels.
%
%
%   Results
%   -------
%
%   Figure 1: Results
%
%   
%
%   References:
%     N. Komodakis and J.-C. Pesquet. Playing with duality: An overview of
%     recent primal-dual approaches for solving large-scale optimization
%     problems. arXiv preprint arXiv:1406.5429, 2014.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/demos/demo_fbb_primal_dual.php

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


%% Initialization

clear;
close all;

init_unlocbox;
ltfatstart;
%%
sigma = 0.1;
missing_ratio = 0.4;

lambda = 0.05;
tau = 0.1;
verbose = 2;

%% Creating the problem

img = barbara();

[nx,ny] = size(img);
A = rand(nx,ny)>missing_ratio;
noisy_img = img + sigma * randn(nx,ny);
b = A .* noisy_img;


%% Defining proximal operators

% setting the function f2 
ffid.grad = @(x) 2 * A .* (A.*x - b);
ffid.eval = @(x) norm(A(:).*x(:)-b(:))^2;
ffid.beta = 2;

% setting the function f1 (norm TV)
param_tv.verbose=verbose - 1;
param_tv.maxit=50;

ftv.prox=@(x, T) prox_tv(x, T*lambda, param_tv);
ftv.eval=@(x) lambda * norm_tv(x);   

% wavelet
Nlevel = 5;
W =@(x) fwt2(x,'db8',Nlevel);
Wt = @(x) ifwt2(x,'db8',Nlevel);
paraml1.verbose = verbose -1;
fw.prox = @(x,T) prox_l1(x,tau*T,paraml1);
fw.eval = @(x) tau*norm(W(x),1);
fw.L =W;
fw.Lt = Wt;

%% Solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 50;        % maximum iteration
%param_solver.gamma = 0.5;       % stepsize (beta is equal to 2)
param_solver.tol = 1e-6;        % Tolerance to stop iterating
fig=figure(100);
param_solver.do_sol=@(x) plot_image(x,fig);
% solving the problem
sol = fb_based_primal_dual(b,ffid,ftv, fw,param_solver);
close(100);

%% displaying the result
imagesc_gray(img, 1, '(a) Original image',221,[0 1]);
imagesc_gray(noisy_img, 1, '(b) Noisy image',222,[0 1]);
imagesc_gray(b, 1, '(c) Measurements',223,[0 1]);
imagesc_gray(sol, 1, '(d) Solution of optimization',224,[0 1]);  
paramplot.position = [100 100 500 500];
gsp_plotfig('inpainting_fbb_primal_dual',paramplot);

%% Closing the toolbox
close_unlocbox();







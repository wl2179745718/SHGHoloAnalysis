function [sol,info] = prox_l0(x, gamma, param)
%PROX_L0 Proximal operator with L0 norm
%   Usage:  sol=prox_l0(x)
%           sol=prox_l0(x, gamma)
%           sol=prox_l0(x, gamma, param)
%           [sol, info]=prox_l0(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   PROX_L0(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * || z ||_0
%
%   param is a Matlab structure containing the following fields:
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.k : number of non zero element (if not defined, it uses gamma
%     to dertermine how many coefficients are kept
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used 
%
%   See also:  proj_b1 prox_l1
%
%   Url: http://unlocbox.sourceforge.net/doc/prox/prox_l0.php

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
% Date: Nov 2012
% Testing: test_prox_l1

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end

% test the parameters
if test_gamma(gamma)
    sol = x;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end




if isfield(param,'k') % k defined
    [~,ind] = sort(x,'ascend');
    sol = x;
    sol(ind(1:(end-param.k))) = 0;
    norm_l0 = param.k;    
else % k not defined
    hard = @(x,T) x .* (abs(x) > T);
    sol = hard(x,gamma);
    norm_l0 = sum(abs(sol(:))>0);
end

crit = 'TOL_EPS';
iter = 1;

% Log after the prox l1
if param.verbose >= 1
    fprintf(['  prox_L1: ||x||_0 = %e,', ...
        ' %s, iter = %i\n'], norm_l0, crit, iter);
end


info.algo=mfilename;
info.iter=iter;
info.final_eval=norm_l0;
info.crit=crit;
info.time=toc(t1);
end




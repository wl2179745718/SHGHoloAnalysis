function [sol, info,objective] = fbf_primal_dual(x_0,f1, f2, f3, param)
%FBF_PRIMAL_DUAL forward backward forward primal dual
%   Usage: sol = fbf_primal_dual(x_0,f1, f2, f3, param);
%          sol = fbf_primal_dual(x_0,f1, f2, f3);
%          [sol,info,objective] = fbf_primal_dual(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f1    : First function to minimize
%         f2    : Second function to minimize
%         f3    : Third function to minimize
%         param : Optional parameters
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objectiv function each iteration)
%
%   FBF_PRIMAL_DUAL (using forward backward forward based primal dual)
%   solves:
%
%      sol = argmin f1(x) + f2(Lx) + f3(x)
%
%   where  x is the optimization variable with f_1 or f_3 a smooth
%   function and L a linear operator. f_1 and f_3 are defined like
%   other traditional functions.
%
%   Note that f2 is a structure of a functions with:
%
%    f2.eval(x_i) : an operator to evaluate the function
%    f2.prox(x_i, gamma) : an operator to evaluate the prox of the function
%
%   Optionally you can define
%
%    f2.L  : linear operator, matrix or operator (default identity)
%    f2.Lt : adjoint of linear operator, matrix or operator (default identity)
%
%   param a Matlab structure containing solver paremeters. See the
%   function SOLVEP for more information. Additionally it contains those
%   aditional fields:  
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         max_i ||  y_i(t) - y_i(t-1) ||  / ||y(t) ||< tol,
%      
%     where  y_i(t) are the dual variable of function i at itertion t*
%     by default, tol=10e-4.
%
%       Warning! This stopping criterion is different from other solvers!
%
%    param.norm_L : bound on the norm of the operator L (default: 1), i.e.
%
%        ` ||L x||^2 <= norm_L^2 * ||x||^2 
%
%   See also: solvep fb_based_primal_dual
%
%   Demos:  demo_fbb_primal_dual
%
%   References:
%     N. Komodakis and J.-C. Pesquet. Playing with duality: An overview of
%     recent primal-dual approaches for solving large-scale optimization
%     problems. arXiv preprint arXiv:1406.5429, 2014.
%     
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/fbf_primal_dual.php

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
 
% Author: Vassilis Kalofolias, Nathanael Perraudin
% Date: June 2015
% Testing: test_solvers

param.algo = 'FBF_PRIMAL_DUAL';
[sol, info,objective] = solvep(x_0,{f1, f2, f3},param);

end

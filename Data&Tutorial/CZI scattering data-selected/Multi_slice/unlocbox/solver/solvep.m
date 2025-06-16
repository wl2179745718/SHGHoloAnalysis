function [sol, info,objective] = solvep(x_0, F, param)
%SOLVEP solve a minimization problem
%   Usage: sol = solvep(x_0, F, param);
%          sol = solvep(x_0, F);
%          [sol,infos,objectiv] = solvep(...);
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         F     : array of function to minimize (structure)
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%         objective: vector (evaluation of the objectiv function each iteration)
%
%   SOLVEP solves:
%
%      sol = argmin sum_i fi(x)     for x belong to R^N
%
%   where x is the variable.
%
%   x_0 is the starting point of the algorithm. A good starting point
%   could significantly reduce the computation time
%
%   F is an array of structure representing convex function to be
%   minimized. Those Functions
%   can be minimized thanks to their gradient (if they are differentiable)
%   or thanks to their proximal operator. As a result the algorithm will needs some of those operators. The easiest way to define a
%   function f1 is to create a struture with the fields f1.eval, f1.grad
%   and f1.prox. Those field all contatains an inline function that
%   compute respectively the evaluation of the function itself, the
%   gradient or the proximal operator. Depending on the solver, not all
%   this operators are necessary. See each solver documentation for
%   details. When three functions are defined, F = {f1, f2, f3}.
%   
%   param a Matlab structure containing the following fields:
%
%    param.gamma : is the step size. Watch out, this parameter is
%     bounded. It should be below 1/beta (*f2 is beta Lipchitz
%     continuous). By default, it is computed with the lipschitz constant
%     of all smooth functions.
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%           (  n(t) - n(t-1) )  / n(t) < tol,
%      
%       where  n(t) = f_1(x)+f_2(x) is the objective function at iteration t*
%       by default, tol=10e-4.
%
%    param.algo : solver used for the problem. Determined
%     automatically with the functions in f.
%
%    param.maxit : is the maximum number of iteration. By default, it is 200.
% 
%    param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%    param.abs_tol : If activated, this stopping critterion is the
%     objectiv function smaller than param.tol. By default 0.
%
%    param.use_dual : If activated, use the norm of the dual variable
%     instead of the evaluation of the function itself for stopping
%     criterion. This is used in ADMM and SDMM for instance. To use it, the
%     algorithm needs to store the dual variable in s.dual_var.
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the objectivs functions
%
%    info.crit : Stopping critterion used 
%
%    info.rel_norm : Relative norm at convergence 
%
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/solvep.php

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
% Date: 22 Feb 2015
% Testing: test_solver

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end

if nargin<2
    error('Not enough input arguments')
end


if ~isfield(param, 'tol'), param.tol=10e-4 ; end
if ~isfield(param, 'maxit'), param.maxit=200; end
if ~isfield(param, 'verbose'), param.verbose=1 ; end
if ~isfield(param, 'lambda'), param.lambda=1 ; end
if ~isfield(param, 'fast_eval'), param.fast_eval = 0  ; end
if ~isfield(param, 'abs_tol'), param.abs_tol = 0  ; end
if ~isfield(param, 'use_dual'), param.use_dual = 0  ; end


% test the input for eval
if ~iscell(F)
    F = {F};
end
F = test_eval(F);


% test input for grad and prox
[Fg,Fp] = prepare_function(F,param);
if ~isfield(param, 'gamma'), 
    if numel(Fg)
        param.gamma = compute_gamma(Fg);
    else
        param.gamma = 1; 
    end
else
    if param.verbose >= 1
        fprintf('The time step is set manually to : %g\n', param.gamma);
    end    
end

if ~isfield(param, 'algo'), param.algo = select_solver(Fg,Fp)  ; end


algo = get_algo(param.algo);

% Transform all smooth functions into one function.
fg = add_smooth_functions(Fg);

if param.verbose>=1, 
    fprintf(['Algorithm selected:', algo.name,' \n']);
end

[sol,s,param] = algo.initialize(x_0, fg, Fp, param);


% Initialization
if param.use_dual && isfield(s,'dual_var')
    [curr_eval, dual_var_old] = eval_dual_var(s.dual_var);
else
    curr_eval = eval_function(fg,Fp,x_0);
end
[~,~,prev_eval,iter,objective,~] = convergence_test(curr_eval);


% Main loop
while 1

    if param.verbose >= 2
        fprintf('Iteration %i:\n', iter);
    end
    [sol, s] = algo.algorithm(x_0, fg, Fp, sol, s, param);
    
    % Global stopping criterion
    if param.use_dual && isfield(s,'dual_var')
        [curr_eval, dual_var_old] = eval_dual_var(s.dual_var,dual_var_old);
    else
        curr_eval = eval_function(fg,Fp,sol,s,param);
    end   
    [stop,rel_eval,prev_eval,iter,objective,crit] = convergence_test(curr_eval,prev_eval,iter,objective,param,s);
    [sol, param] = post_process(sol, iter, curr_eval, prev_eval, objective, param);
    if param.verbose >= 2
        if param.use_dual && isfield(s,'dual_var')
            fprintf('   Dual relative norm = %e\n', curr_eval);
        else
            fprintf('  ||f|| = %e, rel_norm = %e\n', ...
            curr_eval, rel_eval);
        end
    end
    if stop
        break;
    end
    
end

sol = algo.finalize(x_0, fg, Fp, sol, s, param);

% Log
if param.verbose>=2
    % Print norm
    fprintf(['\n ',algo.name,':\n']);
    if param.use_dual && isfield(s,'dual_var')
        fprintf(' Final dual relative evaluation: %e\n', curr_eval );    
    else
        fprintf(' Final relative evaluation: %e\n', rel_eval );   
        fprintf(' ||f|| = %e\n', curr_eval );  
    end
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit);
elseif param.verbose>=1
    if param.use_dual && isfield(s,'dual_var')
        fprintf([algo.name,': Final dual relative evaluation = %e, it = %i, %s\n'], ...
                        curr_eval, iter,crit);
    else
        fprintf([algo.name,': f(x^*) = %e, rel_eval = %e, it = %i, %s\n'], ...
                        curr_eval, rel_eval, iter,crit);
    end
end

% Here if use_dual && isfield(s,'dual_var'), we might return something
% else. To be fixed...

info.algo=algo.name;
info.iter=iter;
info.final_eval=curr_eval;
info.crit=crit;
info.time=toc(t1);
info.rel_norm=rel_eval;

end


function solver = select_solver(Fg,Fp)
    n = numberofL(Fp);
    if numel(Fg)
        if n>2
            error('Sorry, no solver is able to solve your problem yet!')
        end
        if numel(Fp)==0
            solver = 'GRADIENT_DESCENT';
        elseif numel(Fp)==2
            solver = 'FB_BASED_PRIMAL_DUAL';
        elseif (numel(Fp)<=2) && n
            solver = 'FB_BASED_PRIMAL_DUAL';            
        elseif numel(Fp)==1
            solver = 'FORWARD_BACKWARD';
        else
            solver = 'GENERALIZED_FORWARD_BACKWARD';
        end
    else
        if numel(Fp) == 1
            error('Do you really want to minimize only one function?')        
        elseif (numel(Fp)<=2) && (n == 1)
            solver = 'FB_BASED_PRIMAL_DUAL';    
        elseif (n>1)
            solver = 'SDMM';
        elseif numel(Fp) == 2
            solver = 'DOUGLAS_RACHFORD';
        else
            solver = 'PPXA';
        end
    end
        
end

function gamma = compute_gamma(Fg)
    beta = 0;
    for ii = 1:length(Fg)
        beta = beta + Fg{ii}.beta;
    end
    if beta == 0;
        gamma = 1;
    elseif beta >0
        gamma = 1/beta;
    else
        error('Error in the libschitz constant!')
    end
end

function algo = get_algo(name)
    if isstruct(name)
        algo = name;
    elseif ischar(name)
        algo = algoname(name);
    else
        error('The algorithm is not a string and not a struct!')
    end
end


function [rel_norm, x_old] =  eval_dual_var(x1,x2)

    if nargin<2
        if iscell(x1)
            x2 = cell(length(x1),1);
            for ii = 1:length(x1)
                x2{ii} = 0;
            end
        else
            x2 = 0;
        end
    end
    
    if iscell(x1)
        rel_norm = 0;
        for ii = 1:length(x1)
            tmp = norm(x1{ii}(:)-x2{ii}(:))/(norm(x1{ii}(:))+eps);
            if tmp>rel_norm
                rel_norm = tmp;
            end
        end
    else
        rel_norm = norm(x1(:)-x2(:))/(norm(x1(:))+eps);
    end


    x_old = x1;
end

function n = numberofL(Fp)
% Return the number of functions with a linear opeartor inside
n = 0;
    for ii = 1:length(Fp)
        if isfield(Fp{ii},'L') 
            n = n + 1;
            if ~isfield(Fp{ii},'Lt')
                warning('You did not define the Lt operator!');
            end
        end
    end
end

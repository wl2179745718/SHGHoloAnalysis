function s = demo_forward_backward_alg()
%DEMO_FORWARD_BACKWARD_ALG Demonstration to define a personal solver
%   Usage : param.algo = demo_forward_backward_alg();
%
%   This function returns a structure containing the algorithm. You can
%   lauch your personal algorithm with the following:
%
%           param.algo = demo_forward_backward_alg();
%           sol = solvep(x0, {f1, f2}, param);
%
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/demo_forward_backward_alg.php

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

    % This function returns a structure with 4 fields:
    % 1) The name of the solver. This is used to select the solvers.
    s.name = 'DEMO_FORWARD_BACKWARD';
    % 2) A method to initialize the solver (called at the beginning)
    s.initialize = @(x_0, fg, Fp, param) ...
      forward_backward_initialize(x_0,fg,Fp,param);
    % 3) The algorithm itself (called at each iterations)
    s.algorithm = @(x_0, fg, Fp, sol, s, param) ...
      forward_backward_algorithm(fg, Fp, sol, s, param);
    % 4) Post process method (called at the end)
    s.finalize = @(x_0, fg, Fp, sol, s, param) sol;
    % The variables here are
    %   x_0 : The starting point
    %   fg  : A single smooth function 
    %         (if fg.beta == 0, no smooth function is specified) 
    %   Fp  : The non smooth functions (a cell array of structure)
    %   param: The structure of optional parameter
    %   s   : Intern variables or the algorithm
    %   sol : Current solution
end

function [sol, s, param] = forward_backward_initialize(x_0,fg,Fp,param)

    % Handle optional parameter. Here we need a variable lambda.
    if ~isfield(param, 'lambda'), param.lambda=1 ; end

    % All intern variables are stored into the structure s
    s = struct;
    % *sol* is set to the initial points
    sol = x_0;
    
    if numel(Fp)>1
        error(['This solver can not be used to optimize',...
          ' more than one non smooth function']);
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
    
end

function [sol, s] = forward_backward_algorithm(fg, Fp, sol, s, param)
    % The forward backward algorithm is done in two steps
    %  1) x_n = prox_{f, gamma} ( sol - gamma grad_fg(sol) )
    s.x_n = Fp{1}.prox( sol - param.gamma*fg.grad(sol), param.gamma);
    %  2) Updates
    %     sol = sol + lambda * (x_n -sol)
    sol = sol + param.lambda * (s.x_n - sol);
end

function s = forward_backward_alg()
   s.name = 'FORWARD_BACKWARD';
   s.initialize = @(x_0, fg, Fp, param) forward_backward_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) forward_backward_algorithm(fg, Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = forward_backward_initialize(x_0,fg,Fp,param)

    if ~isfield(param, 'method'), param.method='FISTA' ; end
    if ~isfield(param, 'lambda'), param.lambda=1 ; end
    s.method = param.method;
    s.lambda = param.lambda;
    
    s.x_n = {};
    if strcmp(s.method, 'FISTA')
        s.u_n = x_0;
        s.tn = 1;
    end
    sol = x_0;
    if numel(Fp)>1
        error('This solver can not be used to optimize more than one non smooth function')
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
end


function [sol, s] = forward_backward_algorithm(fg, Fp, sol, s, param)
    
    if strcmp(s.method, 'FISTA')

%         % FISTA algorithm
%         s.x_n{1} = Fp{1}.prox_ad(s.u_n-param.gamma*fg.grad(s.u_n), param.gamma);
%         tn1 = (1 + sqrt(1+4*s.tn^2))/2;
%         s.u_n = s.x_n{1}{1} + (s.tn-1)/tn1*(s.x_n{1}{1}-sol);
%         % updates
%         sol = s.x_n{1}{1};
%         s.tn = tn1;
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/alg/forward_backward_alg.php

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
        % FISTA algorithm
        s.x_n{1} = Fp{1}.prox(s.u_n-param.gamma*fg.grad(s.u_n), param.gamma);
        tn1 = (1 + sqrt(1+4*s.tn^2))/2;
        s.u_n = s.x_n{1} + (s.tn-1)/tn1*(s.x_n{1}-sol);
        % updates
        sol = s.x_n{1};
        s.tn = tn1;
    else
%         s.x_n{1} = Fp{1}.prox_ad( sol - param.gamma*fg.grad(sol), param.gamma);
%         sol = sol + param.lambda*(s.x_n{1}{1} - sol);        
        s.x_n{1} = Fp{1}.prox( sol - param.gamma*fg.grad(sol), param.gamma);
        sol = sol + param.lambda*(s.x_n{1} - sol);
    end

end


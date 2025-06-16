function s = chambolle_pock_alg()
   s.name = 'CHAMBOLLE_POCK';
   s.initialize = @(x_0, fg, Fp, param) chambolle_pock_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) chambolle_pock_algorithm(Fp, sol,s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = chambolle_pock_initialize(x_0,fg,Fp,param)
    
    error('This solver contains a bug and is not working yet!')
    
    if ~isfield(param, 'tau'), param.tau=1 ; end
    if ~isfield(param, 'rho'), param.rho=1 ; end
    if ~isfield(param, 'L'), param.L=@(x) x; end
    if ~isfield(param, 'Lt'), param.Lt=@(x) x; end

    s.tau = param.tau;
    s.rho = param.rho;
    if isa(param.L,'numeric')
       s.OpL= @(x) param.L*x;
    else
       s.OpL= param.L;
    end

    if isa(param.Lt,'numeric')
       s.OpLt= @(x) param.Lt*x;
    else
       s.OpLt= param.Lt;
    end



    s.dual_var = s.OpL(x_0);
    sol = x_0;
    s.x_n = x_0;


    if fg.beta
        error('CHAMBOLLE POCK needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('CHAMBOLLE POCK needs exactly 2 functions')
    end
    
    param.abs_tol = 1;
    param.use_dual = 1;
   
    
end


function [sol, s] = chambolle_pock_algorithm(Fp, sol, s, param)
    
    % Algorithm
    s.dual_var = prox_adjoint( s.dual_var + s.rho *s.OpL(sol), s.rho, Fp{2});
    x_n_old = s.x_n;
    s.x_n = Fp{1}.prox( s.x_n + s.tau * s.OpLt(s.dual_var), s.tau);
    sol = s.x_n + param.gamma*(s.x_n - x_n_old);
  
end

%
%   Url: http://unlocbox.sourceforge.net/doc/solver/alg/chambolle_pock_alg.php

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


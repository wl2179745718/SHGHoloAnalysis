function s = sdmm_alg()
   s.name = 'SDMM';
   s.initialize = @(x_0, fg, Fp, param) sdmm_initialize(fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) sdmm_algorithm(Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s,param] = sdmm_initialize(fg,Fp,param)


    if ~isfield(param, 'Qinv'), 
        param.Qinv=@(x) x./numel(Fp); 
        warning('The parameter Qinv is set automatically... It might be wrongly set')
    end

    if isa(param.Qinv,'numeric')
        s.QinvOp= @(x) param.Qinv*x;
    else
        s.QinvOp= param.Qinv;
    end
    s.L = cell(length(Fp),1);
    s.Lt = cell(length(Fp),1);
    s.x_n = cell(length(Fp),1);
    s.u_n = cell(length(Fp),1);
    for ii = 1:length(Fp)
        if ~isfield(Fp{ii}, 'L'), Fp{ii}.L = @(x) x; end
        if ~isfield(Fp{ii}, 'Lt'), Fp{ii}.Lt = @(x) x; end
%         if size(F{i}.L,1)==0, F{i}.L=eye(length(F{i}.x0)) ; end
%         if size(F{i}.x0,2)>size(F{i}.x0,1), F{i}.x0=F{i}.x0'; end
%        s.x_n{ii}{1} = Fp{ii}.x0;
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/alg/sdmm_alg.php

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
        s.x_n{ii} = Fp{ii}.x0;
        s.dual_var{ii} = zeros(size(Fp{ii}.x0));
        
        % Check how the operator is given
        if isa(Fp{ii}.L,'numeric')
            s.OpL{ii}= @(x) Fp{ii}.L*x;
        else
            s.OpL{ii}= Fp{ii}.L;
        end

        if isa(Fp{ii}.Lt,'numeric')
            s.OpLt{ii}= @(x) Fp{ii}.Lt*x;
        else
            s.OpLt{ii}= Fp{ii}.Lt;
        end
    end

    % Initialization
    sol = 0;
    for ii = 1:length(Fp)
%        sol = sol +s.OpLt{ii}( s.x_n{ii}{1} );
        sol = sol +s.OpLt{ii}( s.x_n{ii} );
    end
    
    sol = s.QinvOp(sol);     
    
    if fg.beta
        error('SDMM needs only function with proximal operators')
    end
    
    param.abs_tol = 1;
    param.use_dual = 1;
end


function [sol, s] = sdmm_algorithm(Fp, sol, s, param)    

    for ii = 1 : length(Fp)
        s_n = s.OpL{ii}(sol);
%        s.x_n{ii} = Fp{ii}.prox_ad(s_n + s.dual_var{ii},param.gamma);
       s.x_n{ii} = Fp{ii}.prox(s_n + s.dual_var{ii},param.gamma);
%        s.dual_var{ii} = s.dual_var{ii} + s_n - s.x_n{ii}{1} ;        
        s.dual_var{ii} = s.dual_var{ii} + s_n - s.x_n{ii} ;
    end
    
    sol = 0; 
    for ii = 1 : length(Fp)
%        sol = sol + s.OpLt{ii}(s.x_n{ii}{1} - s.dual_var{ii});
        sol = sol + s.OpLt{ii}(s.x_n{ii} - s.dual_var{ii});
    end
    sol = s.QinvOp(sol);
       
end


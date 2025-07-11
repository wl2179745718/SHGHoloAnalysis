function [Fg,Fp] = prepare_function(F,param)
%PREPARE_FUNCTION test if function do possess correct fields and sort them
%
%
%   Url: http://unlocbox.sourceforge.net/doc/solver/misc/prepare_function.php

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

% Test if the algorithm requires only non smooth functions
smooth = 1;
if isfield(param, 'algo')
    name = lower(getname(param.algo));
    if strcmp(name,'douglas_rachford') || strcmp(name,'ppxa') || ...
            strcmp(name,'admm') || strcmp(name,'sdmm') || ...
            strcmp(name,'chambolle_pock')
        smooth = 0;
    end
end

if ~iscell(F)
    F = {F};
end
    
% number of function
m = length(F);
Fg = {};
Fp = {};
for ii = 1:m
    if isfield(F{ii},'grad') && smooth
        if ~isfield(F{ii},'beta')
            F{ii}.beta = 1;
            warning('Please specify a lipshitz constant of the gradient. 1 will be used instead')
        end
        Fg{length(Fg)+1,1} = F{ii};         %#ok<AGROW>
    elseif isfield(F{ii},'proxL')
        F{ii}.prox_adL = @(x,T) prox_adL(x,T,F{ii},param);
        Fp{length(Fp)+1,1} = F{ii};         %#ok<AGROW>
        if isfield(F{ii},'prox')
            F{ii}.prox_ad = @(x,T) prox_ad(x,T,F{ii},param);
        end
    elseif isfield(F{ii},'prox')
        F{ii}.prox_ad = @(x,T) prox_ad(x,T,F{ii},param);
        Fp{length(Fp)+1,1} = F{ii};         %#ok<AGROW>
    else
        if isfield(F{ii},'grad')
            error('This solver require only smooth function')
        else
            error('No grad or prox function defined');
        end
    end
    
    if isfield(F{ii},'prox') && isfield(F{ii},'grad') && smooth && param.verbose
        fprintf('A function has both the prox and a grad fields. The gradient is used\n');
    end
end


end

function s = prox_ad(x,T,f,param)

if param.fast_eval
    s = {};
    [s{1:2}] = f.prox(x,T);
else
    s = {};
    [s{1}] = f.prox(x,T);
end

end

function s = prox_adL(x,T,f,param)

if param.fast_eval
    s = {};
    [s{1:2}] = f.proxL(x,T);
else
    s = {};
    [s{1}] = f.proxL(x,T);
end

end

function name = getname(algo)
    if ischar(algo)
        name = algo;
    elseif isstruct(algo)
        name = algo.name;
    else
        error('Algo is not a struct or a string')
    end
end



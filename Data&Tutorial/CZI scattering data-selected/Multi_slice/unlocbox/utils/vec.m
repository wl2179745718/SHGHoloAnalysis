function r=vec(x)
%VEC vectorize x
%   Usage: r = vec(x);
%
%   Inputs parameters:
%       x   : vector or matrix
%   Outputs parameters:
%       r   : row vector
%
%   This function vectorize x.
%
%
%   Url: http://unlocbox.sourceforge.net/doc/utils/vec.php

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


r=x(:);


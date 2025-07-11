function [ Lx ] = laplacianx_op( I )
%LAPLACIANX_OP dimentional Laplacian
%   Usage:  [Lx] = laplacianx_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         Lx    : Laplacian along x
%
%   Compute the sum of the laplacian along x. This operator is 
%   self-adjoint.
%
%           Lx = I_xx
%
%   See also: laplacian_op laplaciany_op div_op gradient_op
%
%
%   Url: http://unlocbox.sourceforge.net/doc/utils/laplacianx_op.php

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

% Author: Nathnaael Perraudin
% Date  : 13 September 2013


dx = [I(2:end, :,:)-I(1:end-1, :,:) ; zeros(1, size(I, 2),size(I, 3))];


Lx = [dx(1, :,:) ; dx(2:end-1, :,:)-dx(1:end-2, :,:) ; -dx(end-1, :,:)];


end


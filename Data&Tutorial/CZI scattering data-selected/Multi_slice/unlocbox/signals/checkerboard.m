function [im] = checkerboard()
%CHECKERBOARD Load the 'checkerboard' test signal
%   Usage: im = checkerboard();
%
%   Input parameters:
%       none
%   Output parameters:
%       ime    : image
%
%   Example
%   -------
%   
%   Load the image and display it:
%
%       im = checkerboard();
%       imagescgray(im);
%
%
%   Url: http://unlocbox.sourceforge.net/doc/signals/checkerboard.php

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
% Date: 6 February 2015
 

f = mfilename('fullpath');

% Load the signal

im = imread([f, '.png']);

im = rgb2gray(im);

im = double(im) / 255;


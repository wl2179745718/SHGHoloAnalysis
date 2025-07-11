function [U,S,V] = svdsecon(X,k)
%SVDECON Fast svds when n<<m
%
%   Url: http://unlocbox.sourceforge.net/doc/utils/svdsecon.php

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

    [m,n] = size(X);

    if  m <= n
        X2 = X*X';
        [U,E] = eigs(X2,k);

        [e,ix] = sort(abs(diag(E)),'descend');
        U = U(:,ix);    

        V = X'*U;
        s = sqrt(e);
        
        V = bsxfun(@times, V, 1./s');
        S = diag(s);
    else
        X2 = X'*X; 
        [V,E] = eigs(X2,k);

        [e,ix] = sort(abs(diag(E)),'descend');
        V = V(:,ix);    

        U = X*V; 
        s = sqrt(e);
        U = bsxfun(@times, U, 1./s');
        S = diag(s);
    end
end

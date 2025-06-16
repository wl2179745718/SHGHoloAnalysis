function J = ric_besselj_derivative(nu, x, flag)
% J = ric_besselj_derivative(nu, x, flag) use the recursive
% relationship to calculate the first derivative of the Riccati-Bessel
% funtion.
%
% The Riccati-Bessel function is defined as
%
% J_{nu}(x) = \sqrt{\frac{\pi x}{2}} J_{nu+1/2}(x)
%
% Input:
%
% nu         order of the Riccati-Bessel function.  Must be a
%            column vector.
%
% x          Must be a row vector.
%
% flat       1, first order derivative order; 2, second order derivative
%

    if (nargin == 2)
        flag = 1;
    end

    x    = reshape(x, 1, length(x));
    nu   = reshape(nu, length(nu) ,1);

    temp = ones(length(nu), 1)*x;
    if (flag == 1)
        J = 0.5*(           ric_besselj(nu-1, x) ...
                 + 1./temp.*ric_besselj(nu,   x) ...
                 -          ric_besselj(nu+1, x));
    elseif (flag ==2)
        J = 0.5*(              ric_besselj_derivative(nu-1, x) ...
                 +    1./temp.*ric_besselj_derivative(nu,   x) ...
                 - temp.^(-2).*ric_besselj(nu, x) ...
                 -             ric_besselj_derivative(nu+1, x));
    else
        error('This script only handles first and second derivative.')
    end
    

    temp2 = ones(length(nu),1)*x;
    J(temp2==0) = 0; % x = 0, all zeros
    
    temp1 = nu*ones(1, length(x));
    if (flag ==1)
        J(temp1==0 & temp2==0) = 1;
    end
end


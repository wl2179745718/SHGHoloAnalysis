function J = ric_besselj(nu,x)
% J = ric_besselj(nu, x) implement the Riccati-Bessel functions.
% 
% J_{nu}(x) = \sqrt{\frac{\pi x}{2}} J_{nu+1/2}(x)
% 
% Input:
%
% nu  order of the Bessel function. Must be a column vector.
% x   a row of complex vectors
%

    nu = reshape(nu, length(nu), 1); % impose to be a column vector
    x  = reshape(x, 1, length(x)); % impose x to be a row vector
    
    a = zeros(length(nu), length(x));
    for iNu = 1:length(nu)
        a(iNu, :) = besselj(nu(iNu)+0.5,x);
    end
    
    J = sqrt(pi/2.*(ones(length(nu), 1)*x)).*a;



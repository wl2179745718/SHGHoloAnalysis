function Y = ric_bessely(nu,x)
% Y = ric_bessely(nu, x) implements the Riccati-Neumann's functions
%
% Y_{nu}(x) = \sqrt{\frac{\pi x}{2}} Y_{nu+1/2}(x)
%
% where
% 
% nu  order of the Bessel's function. Must be a column vector.
%
% x   must be a row complex vector.
%


    nu = reshape(nu, length(nu) ,1);
    x  = reshape(x, 1, length(x));

    a = zeros(length(nu), length(x));
    for iNu = 1:length(nu)
        a(iNu, :) = bessely(nu(iNu)+0.5,x);
    end

    Y = sqrt(pi/2.*(ones(length(nu), 1)*x)).*a;
    
    if (sum(find(x==0))~=0)
        disp('ric_bessely evaluated at x=0. Return -inf');
    end
    
    temp2 = ones(length(nu),1)*x;
    Y(temp2==0) = -inf;
    
    temp1 = nu*ones(1, length(x));
    Y(temp1==0 & temp2==0) = -1;
    
    
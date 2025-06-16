function H = besselh_derivative(nu, K, x)
% H = besselh_derivative(nu, K, x) use the recursive relationship to
% calculate the derivative of the bessel funtion of the third kind,
% also known as the Hankel function.
%
% nu      order of the bessel function
%
% K = 1   if it is Hankel function of the first kind; K=2 if it is
%         Hankel function of the second kind
%
% x       real or complex scalar
%
%
    
    if (K ~=1 && K~=2)
        error('Improper kind of Hankel function');
    end
    
    H = 0.5*(besselh(nu-1,K,x)-besselh(nu+1,K,x));
end 



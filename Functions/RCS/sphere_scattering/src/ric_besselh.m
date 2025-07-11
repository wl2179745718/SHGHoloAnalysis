function H = ric_besselh(nu,K,x)
% H = ric_besselh(nu, K, x) implement the Riccati-Bessel function,
% which is defined as
%
% H_{nu}(x) = \sqrt{\frac{\pi x}{2}} H_{nu+1/2}(x)
%
% Input:
%
% nu  order of the spherical Bessel's function. Must be a column vector.
%
% K   1 for Hankel's function of the first kind; 2 for Hankel's
%     function of the second kind.
%
% x   Must be a row vector.
%

    if (K ~=1 && K ~=2)
        error('Improper kind of Hankel function');
    end
    
    if (K==1)
        H = ric_besselj(nu,x) + 1j*ric_bessely(nu,x);
    else
        H = ric_besselj(nu,x) - 1j*ric_bessely(nu,x);
    end
    

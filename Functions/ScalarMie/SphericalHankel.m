function [h] = SphericalHankel(l,ka)
% spherical Hankel
% https://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html

h = sqrt(pi/2./ka).*besselh(l+1/2,ka);

end
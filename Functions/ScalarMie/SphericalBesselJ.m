function [j] = SphericalBesselJ(l,ka)
% spherical bessel J
% https://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html

j = sqrt(pi/2./ka).*besselj(l+1/2,ka);

end
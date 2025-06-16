function [dh] = dSphericalHankel(l,ka)
% derivative of spherical Hankel

dh = l./ka.*SphericalHankel(l,ka)-SphericalHankel(l+1,ka);

end
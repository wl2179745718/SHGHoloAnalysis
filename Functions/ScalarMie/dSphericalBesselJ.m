function [dj] = dSphericalBesselJ(l,ka)
% derivative of spherical bessel J

dj = l./ka.*SphericalBesselJ(l,ka)-SphericalBesselJ(l+1,ka);

end
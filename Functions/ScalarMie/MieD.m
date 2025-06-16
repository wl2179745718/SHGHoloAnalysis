function [D] = MieD(n0,ns,k0,a,l)
% Coefficient D of Mie scattering

D = ( n0.*SphericalBesselJ(l,ns*k0.*a).*dSphericalBesselJ(l,n0*k0.*a)-ns.*SphericalBesselJ(l,n0*k0.*a).*dSphericalBesselJ(l,ns*k0.*a) )./ ( ns.*SphericalHankel(l,n0*k0.*a).*dSphericalBesselJ(l,ns*k0.*a)-n0.*dSphericalHankel(l,n0*k0.*a).*SphericalBesselJ(l,ns*k0.*a) );

end
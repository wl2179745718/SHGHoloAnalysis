function [U] = ScalarMie(lmax,theta,r,n0,ns,k0,a)
% Scalar Mie

U = zeros([size(r),lmax]);

parfor l=0:(lmax-1)
    l
    Ul(:,:,l+1)=1i^l*(2*l+1).*MieD(n0,ns,k0,a,l).*SphericalHankel(l,n0*k0.*r).*legendreP(l,cos(theta));

end

U = sum(Ul,3);

end
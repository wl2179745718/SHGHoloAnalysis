function [E, inv_phase] = Mie_x_radiation2(n1, n2, d, lambda, Hsize, dpix, z, cen, k_dir, digit)
%   INPUTS
%   n1, refractive index of the immersion medium
%   n2, refractive index of the sphere
%   d, diameter of the particle
%   lambda, wavelength of the light in the ambient medium
%   Hsize, Number of pixels of the recorded image
%   z, propagation distance
%   dpix, pixel size
%   cen: center of the particle (units: pixels)
%   k_dir: illumination direction ([x,y,z], will be normalized inside)
%   digit, precision digit
%
%   OUTPUTS:
%   E = [Ex; Ey; Ez], in Cartesian Coordinates
% last modified by Lei Tian, Dec 15, 2010
% version 3
%   fixed the appoximation in theta. 
% version 4
% 02/27/2013, Lei Tian
% allow off-center particles
%%
% The ONLY approximation is to assume radiating field (kr>>1), which is
% commonly satisfied for visible range.

%%
if nargin == 9
    digit = 6;
end
k_dir = k_dir/norm(k_dir);
pol_dir = [1,0,0] - k_dir(1)*k_dir;
pol_dir = pol_dir/norm(pol_dir);
m = n2/n1;
% wave number
k = 2*pi*n1/lambda;
% size parameter
alpha = 2*pi*n1/lambda*d/2;
% index of refraction
m1=real(m); m2=imag(m);
% coordiates
x = [-ceil(Hsize/2-1)*dpix:dpix:ceil(Hsize/2)*dpix]+cen(1)*dpix;
y = [-ceil(Hsize/2-1)*dpix:dpix:ceil(Hsize/2)*dpix]+cen(2)*dpix;
[xmesh, ymesh] = meshgrid(x,y);
rmesh = sqrt(xmesh.^2+ymesh.^2);

% distance
r = sqrt(rmesh.^2+z^2);
% spherical wave terms
spherical_term = 1i./(k*r).*exp(1i*k*r);
% original scattering angle 
% Note: u = (r_vector.k_dir)/r
u = (xmesh*k_dir(1)+ymesh*k_dir(2)+z*k_dir(3))./r;
% scattering angle
theta = acos(u);
% azimuthal angle
phi = acos((xmesh*pol_dir(1)+ymesh*pol_dir(2)+z*pol_dir(3))./r./sqrt(1-u.^2));
phi(~isfinite(phi)) = 0;
assert(norm(imag(phi),'fro')<1e-4);
phi = real(phi);
% approximate theta, precision up to 5 digits
thetaapprox = thetaround(theta,digit);
[thetacompute, ~, idx] = unique(thetaapprox);
idx = reshape(idx, Hsize, Hsize);
% approximated scattering angle u
ucompute = cos(thetacompute);
num = length(ucompute);

% uapprox = cos(thetaapprox);
% S1 = zeros(Hsize);
% S2 = zeros(Hsize);
% wb=waitbar(0,'Calculating Mie Series...');
S12tmp = zeros(num,2);
parfor j = 1:num
    % S12
    S12tmp(j,:) = Mie_S12(m,alpha,ucompute(j));
    % idx = find(uapprox==ucompute(j));
    % S1(idx) = S12tmp(j,1);
    % S2(idx) = S12tmp(j,2);
%     if mod(j,1000)==0
%         waitbar(j/num,wb);
%     end
end
tmp = S12tmp(:,1);
S1 = tmp(idx);
tmp = S12tmp(:,2);
S2 = tmp(idx);

% close(wb);

% x-polarized light, 
% REF: [1] Ye Pu and Hui Meng, "Intrinsic aberrations due to Mie scattering in
% particle holography," J. Opt. Soc. Am. A 20, 1920-1932 (2003) 
% [2] :/ Tsamg. K. A. Kong, Scattering of Electromagnetic waves
% Ex
E(:,:,1) = spherical_term .* (cos(phi).^2.*u.*S2 +sin(phi).^2.*S1);
% Ey
E(:,:,2) = spherical_term .* (sin(phi).*cos(phi)).*(S2.*u -S1);
% Ez
E(:,:,3) = -spherical_term .* cos(phi).*sin(theta).*S2;

inv_phase = exp(-1i*k*(u.*r));

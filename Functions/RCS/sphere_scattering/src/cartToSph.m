function [r, theta, phi] = cartToSph(x, y, z)
% [r, theta, phi] = cartToSph(x, y, z) converts the cartesian
% coordinate system to the spherical coordinate system according to
% the following definition:
% 
% r        distance from the origin to the point in the interval 
%          [0, \infty)
%
% theta    elevation angle measured between the positive z-axis and
%          the vector in the interval [0, pi]
%
% phi      azimuth angle measured between the positive x-axis and
%          the vector in the interval [0, 2*pi)
%
% 
    
    if (x==0 && y==0 && z==0)
        r     = 0;
        theta = 0;
        phi   = 0;
    else
        hypotxy = hypot(x,y);
        r       = hypot(hypotxy,z);
        theta   = acos(z/r);
        phi     = 0;
        if (x == 0 && y == 0)
            phi = 0;
        elseif (x >= 0 && y >= 0)
            phi = asin(y/hypotxy);
        elseif (x <= 0 && y >= 0)
            phi = pi - asin(y/hypotxy);
        elseif (x <= 0 && y <= 0)
            phi = pi - asin(y/hypotxy);
        elseif (x >= 0 && y <= 0)
            phi = 2*pi + asin(y/hypotxy);
        end
    end
    
    assert(0<=r);
    assert(0<=theta && theta <= pi);
    assert(0<=phi && phi <2*pi);
    
    
function [G] = ScalarG(X, Y, z, k0)
% sclalar Green's function

G = -exp(1i*k0*sqrt(z^2+ X.^2 + Y.^2))./sqrt(z^2+ X.^2 + Y.^2)./pi./4;

end
function [Gk] = GOlivier(X, Y, dx, dy, dz, k0)
% Olivier''s Green function, in k space

%dx = delta(1); dy = delta(2); dz = delta(3); 

G = -exp(1i*k0*sqrt(dz^2+ X.^2 + Y.^2))./sqrt(dz^2+ X.^2 + Y.^2);

Gk = myfft(G)/pi/4;

end
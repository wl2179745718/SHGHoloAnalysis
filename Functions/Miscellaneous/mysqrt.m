function [y] = mysqrt(dksq,z)

y = real(sqrt(dksq)) + 1i*sign(z)*imag(sqrt(dksq));

end
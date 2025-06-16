function [Fx] = myfft(x)

Fx = fftshift(fft2(ifftshift(x)));

end
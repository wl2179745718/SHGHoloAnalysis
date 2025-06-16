function [E] = BPM(E_sca, Z, NA, k, kxx, kyy)

SquareRt=@(a) abs(real(sqrt(a)))+1i*sign(Z)*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field
prop_Z= exp(1i*Z*SquareRt(k^2-(kxx.^2+kyy.^2)));%+1i*eps %Add small absorption term to avoid indeterminates in the angular greens function
%prop_Z(kxx.^2 + kyy.^2 > k^2) = 0;
e = fft2(E_sca).*prop_Z;
e(kxx.^2 + kyy.^2 > (NA*k)^2) = 0;
E = ifft2(e);

end
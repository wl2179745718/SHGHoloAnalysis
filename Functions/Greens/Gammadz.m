function [Gamma] = Gammadz(lambda,fxx,fyy)
% P(kx,ky)
SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field
Gamma=SquareRt((1/lambda)^2-(fxx.^2+fyy.^2));

end
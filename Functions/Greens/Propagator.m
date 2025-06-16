function [angprop] = Propagator(lambda,FX,FY,z)
% P(kx,ky)

%SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field
%prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2));%+1i*eps %Add small absorption term to avoid indeterminates in the angular greens function
%propdz=exp(dz*prop_phs);

angprop = exp(1i*2*pi*z*mysqrt((1/lambda)^2-FX.^2-FY.^2,z));
%propdz = myifft(myfft(e).*angprop);

end
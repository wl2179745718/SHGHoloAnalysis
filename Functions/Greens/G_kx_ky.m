function [AG] = G_kx_ky(fxx,fyy,n_imm,lambda,dz,dGk,Eps)
% The scalar Green's function in the (kx,ky) space

bound = dGk * ( max(max(sqrt(fxx.^2 + fyy.^2))) - n_imm/lambda );

SquareRt=@(a) abs(real(sqrt(a)))+1i*abs(imag(sqrt(a))); % making sure imaginary and real part is pos otherwise will get gain or backward propagating field

%prop_phs= 1i*2*pi*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2));%+1i*eps %Add small absorption term to avoid indeterminates in the angular greens function
%prop=@(z) exp(prop_phs*z);
%Mask=(n_imm/lambda)^2>1.01*((fxx.^2+fyy.^2));

AG= ((-1i.*Propagator(lambda,fxx,fyy,dz)./(4.*pi.*SquareRt((n_imm/lambda)^2-(fxx.^2+fyy.^2)+1i*Eps)))); % Angular Greens function
AG(fxx.^2 + fyy.^2 > (n_imm/lambda + bound)^2) = 0;

end
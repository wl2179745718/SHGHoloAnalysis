function [E] = dyadic_Green(r,k,d)
% Following Eq. 8 of Short-distance expansion for the electromagnetic half-space Green's tensor: general results and an application to radiative lifetime computations
% https://iopscience-iop-org.proxy2.library.illinois.edu/article/10.1088/1751-8113/42/27/275203/pdf
% I believe Vadim has missed a minus sign

d = d';

R = norm(r);
E = exp(1i*k*R)*( ( k^2/R+1i*k/R^2-1/R^3 )*d + ( -k^2/R-3i*k/R^2+3/R^3 )/R^2*tensorprod(r,r',2,1)*d );

end
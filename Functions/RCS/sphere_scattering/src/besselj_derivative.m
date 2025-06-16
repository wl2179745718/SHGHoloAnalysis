function J = besselj_derivative(nu, x)
% J = besselj_derivative(nu,K,x) use the recursive relationship to
% calculate the derivative of the bessel funtion of the first kind.
%
% nu    order of the bessel function
% x     real or complex scalar
%
%
    
    J = 0.5*(besselj(nu-1, x)-besselj(nu+1, x));
end


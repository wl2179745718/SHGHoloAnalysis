function [ alphal ] = Calc_alphal( l )
% Eq. 3, but theta_s has only one value, so no need to integrate. 
% Lang Wang 20230901

alphal = 1i.^l.*(2.*l + 1);

end
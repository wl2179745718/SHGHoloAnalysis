function [Gs] = scalar_Green(k,R)
Gs = exp(1i*k*R)./4./pi./R;
end
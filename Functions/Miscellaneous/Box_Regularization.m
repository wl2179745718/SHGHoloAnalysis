function [Box,N] = Box_Regularization(Box,dr)

N   = round(Box./dr);
N   = N-mod(N,4);
Box = N.*dr;

end
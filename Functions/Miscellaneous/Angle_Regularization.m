function [thetain,fxin] = Angle_Regularization(thetain0,lambda,dfx)

% only define and x-directed spatial frequency
fxin0 = sin(thetain0)/lambda; %input spatial frequency

% set input spatial frequency to the nearest interger multiple of the
% spatial frquency grid for the FFT
fxin = round(fxin0/dfx)*dfx; % actual input spatial frequency
thetain = asin(lambda*fxin); % actual input angle

end
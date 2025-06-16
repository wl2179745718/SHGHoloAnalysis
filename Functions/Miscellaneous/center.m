function [x,y] = center(X,Y,E)
% Calculate the center of a 2D distribution

x = sum(abs(E).*X,'all')./sum(abs(E),'all');
y = sum(abs(E).*Y,'all')./sum(abs(E),'all');

end
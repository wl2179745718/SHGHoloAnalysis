function [theta, R] = XYz2thetaR(d, X, Y, z)

X_1d = reshape(X, [1 size(X,1)*size(X,2)]);
Y_1d = reshape(Y, [1 size(Y,1)*size(Y,2)]);
theta = 0*X_1d;

for ii = 1:size(theta,2)
    V2 = [X_1d(ii) Y_1d(ii) z];
    theta(ii) = atan2(norm(cross(d,V2)),dot(d,V2));
end
R = sqrt(z^2 + X.^2 + Y.^2);

end
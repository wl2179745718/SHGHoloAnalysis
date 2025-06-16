function [theta, R] = XZy2thetaR(d, X, y, Z)

X_1d = reshape(X, [1 size(X,1)*size(X,2)]);
Z_1d = reshape(Z, [1 size(Z,1)*size(Z,2)]);
theta = 0*X_1d;

for ii = 1:size(theta,2)
    V2 = [X_1d(ii) y Z_1d(ii)];
    theta(ii) = atan2(norm(cross(d,V2)),dot(d,V2));
end
R = sqrt(y^2 + X.^2 + Z.^2);

end
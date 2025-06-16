function [U_Mie] = ScalarMieField(X,Y,z,thetain,lmax,n_imm,nsphere,k0,rad)
% the field at the last layer by scalar Mie

X_1d = reshape(X, [1 size(X,1)*size(X,2)]);
Y_1d = reshape(Y, [1 size(Y,1)*size(Y,2)]);
theta1D = 0*X_1d;
V1 = [tan(thetain) 0 1];
for ii = 1:size(theta1D,2)
    V2 = [X_1d(ii) Y_1d(ii) z(end)];
    theta1D(ii) = atan2(norm(cross(V1,V2)),dot(V1,V2));
end
theta_out = reshape(theta1D,size(X));
r_out = sqrt(X.^2 + Y.^2 + z(end)^2);

U_Mie = ScalarMie(lmax,theta_out,r_out,n_imm,nsphere,k0,rad);

end
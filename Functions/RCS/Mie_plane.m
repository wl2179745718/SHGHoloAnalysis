function [E_sca_Mie] = Mie_plane(X,Y,Z,phi,k,a,f,erb,urb,erp,urp,N)
% calculate the field scattered from a sphere to an observation plane

X_1d = reshape(X, [1 size(X,1)*size(X,2)]);
Y_1d = reshape(Y, [1 size(Y,1)*size(Y,2)]);
theta1D = 0*X_1d;
V1 = [tan(phi) 0 1];
for ii = 1:size(theta1D,2)
    V2 = [X_1d(ii) Y_1d(ii) Z];
    theta1D(ii) = atan2(norm(cross(V1,V2)),dot(V1,V2));
end

[~,~,~,ETheta1D] = mieHKURCS(a,f,erb,urb,erp,urp,N,theta1D);
ETheta = reshape(ETheta1D.*cos(theta1D), [size(X,1) size(X,2)]);%.*cos(theta1D)

R = sqrt(X.^2 + Y.^2 + Z^2);
E_sca_Mie = ETheta./( exp(1i*k*Z)/Z ).*( exp(1i*k*R)./R ); %./( exp(1i*k*z(end))/z(end) )
E_sca_Mie = E_sca_Mie./max(max(E_sca_Mie));

end
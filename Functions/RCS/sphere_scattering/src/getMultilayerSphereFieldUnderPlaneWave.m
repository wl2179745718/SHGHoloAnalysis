function [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getMultilayerSphereFieldUnderPlaneWave(radius, ...
                                                      sphere, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
% [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getMultilayerSphereFieldUnderPlaneWave(radius,
%                                              sphere, ...
%                                              background, ....
%                                              sensor_location, ...
%                                              frequency) 
%
%
% Calculate the field as a plane wave is scattered by a
% multi-layer sphere centered at the origin. The incident plane
% wave is polarized in the +x direction propagating in the +z
% direction.
%
%
%    X
%   ^
%   |
%   |
%   |-----> Z
%
%
%   background      k_1        k_2        k_3        k_4        
%   -----------|----------|----------|----------|---------- 0
%             r_1        r_2        r_3        r_4
%
% See Fig. 11-25 in [Balanis1989] for the exact coordinate convention.
%
%
% Input:
%
% radius        Mx1 vector to specify the radii from the outter most
%               layer to the inner most layer. The elements in
%               this variable must be in the descending order. (m)
%
% sphere        Mx1 cell array of objects of DielectricMaterial to
%               specify the dielectric properties from the outter
%               most layer to the inner most layer.
%
% background    Object of DielectricMaterial
%
% sensor_location  3x1 vector (m)
%
% frequency     Nx1 vector (Hz)
%
% threshold     Threshold of the reciprocal condition number as we invert
%               the matrix to calculate the coefficients. Default to eps.
%
%
% Output:
%
% E_r               Nx1 vector (V/m)
%
% E_phi             Nx1 vector (V/m)
%
% E_theta           Nx1 vector (V/m)
%
% H_r               Nx1 vector (A/m)
%
% H_phi             Nx1 vector (A/m)
%
% H_theta           Nx1 vector (A/m)

    p = inputParser;
    p.addRequired('radius', @isvector);
    p.addRequired('sphere', @iscell);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParameter('threshold', eps, @isnumeric);
    p.addParameter('info', 0, @(x)x==0||x==1);
    p.parse(radius, sphere, background, sensor_location, frequency, varargin{:});
    if (p.Results.info)
        disp(p.Results);
    end
    
    MU_0    = 4*pi*1e-7;
    
    nFreq   = length(frequency);
    nSphere = length(sphere);
    
    eta_m   = getIntrinsicImpedance(background, frequency);    
    k_m     = getWaveNumber(background, frequency);
    mu_m    = getComplexPermeability(background, frequency)*MU_0;
        
    
    eta     = zeros(nFreq, nSphere);    
    k_s     = zeros(nFreq, nSphere);
    mu      = zeros(nFreq, nSphere);
    for idx = 1:nSphere
        eta(:,idx) = getIntrinsicImpedance(sphere{idx}, frequency);    
        k_s(:,idx) = getWaveNumber(sphere{idx}, frequency);
        mu (:,idx) = getComplexPermeability(sphere{idx}, frequency)*MU_0;
    end
    
    N       = getNMax(p.Results.radius, ...
                       p.Results.sphere, ...
                       p.Results.background, ...
                       p.Results.frequency);

    % N is a function of frequency. Set to the maximum of N.
    N       = max(N);
    nu      = 1:N;
    
    [r, theta, phi] = cartToSph(sensor_location(1),sensor_location(2),sensor_location(3));

    
    % aux0 denotes the expression
    % assoc_legendre(nu,1,cos(theta))/sin(theta). Here I am using a
    % recursive relation to compute aux0, which avoids the numerical
    % difficulty when theta == 0 or PI.
    aux0    = zeros(1, length(nu)); 
    aux0(1) = -1; 
    aux0(2) = -3*cos(theta);  
    for n = 2:nu(end-1) 
        aux0(n+1) = (2*n+1)/n*cos(theta)*aux0(n) - (n+1)/n*aux0(n-1);
    end
    
    % aux1 denotes the expression
    % sin(theta)*assoc_legendre_derivative(nu,1,cos(theta)).  Here I am
    % also using a recursive relation to compute aux1 from aux0,
    % which avoids numerical difficulty when theta == 0 or PI.
    aux1    = zeros(1, length(nu));
    aux1(1) = cos(theta); 
    for n = 2:nu(end)
        aux1(n) = (n+1)*aux0(n-1)-n*cos(theta)*aux0(n);
    end

    aux1    = ones(nFreq,1)*aux1;
    aux0    = ones(nFreq,1)*aux0;
    
    
    a_n_scat = zeros(nFreq, N);
    a_n      = zeros(nFreq, N, nSphere);
    c_n      = zeros(nFreq, N, nSphere);
    
    b_n_scat = zeros(nFreq, N);
    b_n      = zeros(nFreq, N, nSphere);
    d_n      = zeros(nFreq, N, nSphere);

    I = sqrt(-1);
    for iN           = 1:N;
        x            = I.^(-iN).*(2*iN+1)./(iN.*(iN+1)); % (11-231a) in [Balanis1989]
        
        % Compute coefficient a_n and c_n
        rhs          = zeros(nFreq, nSphere*2); 
        rhs(:,1)     = x*transpose(ric_besselj_derivative(iN,k_m*radius(1)));
        rhs(:,2)     = x*transpose(ric_besselj(iN,k_m*radius(1)));
        
        matrix       = zeros(nFreq, nSphere*2, nSphere*2);
        if (nSphere == 1)
            matrix(:,      1,      1) = -transpose(ric_besselh_derivative(iN,2,k_m*radius));
            matrix(:,      1,      2) =  k_m./k_s.*transpose(ric_besselj_derivative(iN,k_s*radius));
            matrix(:,      2,      1) = -transpose(ric_besselh(iN,2,k_m*radius));
            matrix(:,      2,      2) =  mu_m./mu.*transpose(ric_besselj(iN,k_s*radius));
        else %nSphere >=2
            idx = 1;
            matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_m*radius(idx)));
            matrix(:,idx*2-1,idx*2  ) =  k_m./k_s(:,idx).*transpose(ric_besselh_derivative(iN,1,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2-1,idx*2+1) =  k_m./k_s(:,idx).*transpose(ric_besselh_derivative(iN,2,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh(iN,2,k_m*radius(idx)));
            matrix(:,idx*2  ,idx*2  ) =  mu_m./mu(:,idx).*transpose(ric_besselh(iN,1,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2+1) =  mu_m./mu(:,idx).*transpose(ric_besselh(iN,2,k_s(:,idx)*radius(idx)));
            if (nSphere > 2)
                for idx = 2:nSphere-1
                    matrix(:,idx*2-1,idx*2-2) = -transpose(ric_besselh_derivative(iN,1,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2-1,idx*2  ) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselh_derivative(iN,1,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2-1,idx*2+1) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselh_derivative(iN,2,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2  ,idx*2-2) = -transpose(ric_besselh(iN,1,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh(iN,2,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2  ,idx*2  ) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselh(iN,1,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2  ,idx*2+1) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselh(iN,2,k_s(:,idx)*radius(idx)));
                end
            end
            idx = nSphere;
            matrix(:,idx*2-1,idx*2-2) = -transpose(ric_besselh_derivative(iN,1,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_s(:,idx-1))*radius(idx));
            matrix(:,idx*2-1,idx*2  ) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselj_derivative(iN,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2-2) = -transpose(ric_besselh(iN,1,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh(iN,2,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2  ,idx*2  ) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselj(iN,k_s(:,idx)*radius(idx)));
        end
        
        for iFreq = 1:nFreq
            A                   = squeeze(matrix(iFreq,:,:));
            b                   = transpose(rhs(iFreq,:));
            if (rcond(A) > p.Results.threshold)
                temp                = A\b;
                a_n_scat(iFreq, iN) = temp(1);
                a_n(iFreq, iN, 1:end-1) = temp(3:2:end-1);
                c_n(iFreq, iN, 1:end-1) = temp(2:2:end-2);
                a_n(iFreq, iN, end)     = temp(end);
            end
        end%iFreq
        
        
        % Compute coefficient b_n and d_n
        rhs          = zeros(nFreq, nSphere*2);
        rhs(:,1)     = x*transpose(ric_besselj(iN,k_m*radius(1)));
        rhs(:,2)     = x*transpose(ric_besselj_derivative(iN,k_m*radius(1)));
        
        matrix       = zeros(nFreq, nSphere*2, nSphere*2);
        if (nSphere == 1)
            matrix(:,      1,      1)  = -transpose(ric_besselh(iN,2,k_m*radius));
            matrix(:,      1,      2)  =  k_m./k_s.*transpose(ric_besselj(iN,k_s*radius));
            matrix(:,      2,      1)  = -transpose(ric_besselh_derivative(iN,2,k_m*radius));
            matrix(:,      2,      2)  =  mu_m./mu.*transpose(ric_besselj_derivative(iN,k_s*radius));
        else %nSphere >=2
            idx = 1;
            matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh(iN,2,k_m*radius(idx)));
            matrix(:,idx*2-1,idx*2  ) =  k_m./k_s(:,idx).*transpose(ric_besselh(iN,1,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2-1,idx*2+1) =  k_m./k_s(:,idx).*transpose(ric_besselh(iN,2,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_m*radius(idx)));
            matrix(:,idx*2  ,idx*2  ) =  mu_m./mu(idx).*transpose(ric_besselh_derivative(iN,1,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2+1) =  mu_m./mu(idx).*transpose(ric_besselh_derivative(iN,2,k_s(:,idx)*radius(idx)));
            if (nSphere > 2)
                for idx = 2:nSphere-1
                    matrix(:,idx*2-1,idx*2-2) = -transpose(ric_besselh(iN,1,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh(iN,2,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2-1,idx*2  ) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselh(iN,1,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2-1,idx*2+1) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselh(iN,2,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2  ,idx*2-2) = -transpose(ric_besselh_derivative(iN,1,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_s(:,idx-1)*radius(idx)));
                    matrix(:,idx*2  ,idx*2  ) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselh_derivative(iN,1,k_s(:,idx)*radius(idx)));
                    matrix(:,idx*2  ,idx*2+1) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselh_derivative(iN,2,k_s(:,idx)*radius(idx)));
                end
            end
            idx = nSphere;
            matrix(:,idx*2-1,idx*2-2) = -transpose(ric_besselh(iN,1,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2-1,idx*2-1) = -transpose(ric_besselh(iN,2,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2-1,idx*2  ) =  k_s(:,idx-1)./k_s(:,idx).*transpose(ric_besselj(iN,k_s(:,idx)*radius(idx)));
            matrix(:,idx*2  ,idx*2-2) = -transpose(ric_besselh_derivative(iN,1,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2  ,idx*2-1) = -transpose(ric_besselh_derivative(iN,2,k_s(:,idx-1)*radius(idx)));
            matrix(:,idx*2  ,idx*2  ) =  mu(:,idx-1)./mu(:,idx).*transpose(ric_besselj_derivative(iN,k_s(:,idx)*radius(idx)));
        end
        
        for iFreq = 1:nFreq
            A                   = squeeze(matrix(iFreq,:,:));
            b                   = transpose(rhs(iFreq,:));
            if (rcond(A) > p.Results.threshold)
                temp                = A\b;
                b_n_scat(iFreq, iN) = temp(1);
                b_n(iFreq, iN, 1:end-1) = temp(3:2:end-1);
                d_n(iFreq, iN, 1:end-1) = temp(2:2:end-1);
                b_n(iFreq, iN, end)     = temp(end);
            end
        end%iFreq
    end %iN

    % Due to the large number of orders to evaluate, we may run
    % into numerical difficulty. In this case, we set the
    % cofficients to 0 by default.
    a_n_scat(isnan(a_n_scat)) = 0;
    b_n_scat(isnan(b_n_scat)) = 0;
    a_n(isnan(a_n)) = 0;
    b_n(isnan(b_n)) = 0;
    c_n(isnan(c_n)) = 0;
    d_n(isnan(d_n)) = 0;
    
    if (r <= radius(end)) 
        % Inner most sphere, including the boundary

        x         = k_s(:,end)*r;
        
        alpha     = (transpose(ric_besselj_derivative(nu,x,2))+transpose(ric_besselj(nu,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        alpha(isnan(alpha)) = 0;
        E_r       = -I*cos(phi)*sum(a_n(:,:,end).*alpha, 2);
        H_r       = -I*sin(phi)*sum(b_n(:,:,end).*alpha, 2)./eta(:,end);
        
        alpha     = transpose(ric_besselj_derivative(nu,x)).*aux1;
        beta      = transpose(ric_besselj(nu,x)).*aux0;
        alpha(isnan(alpha)) = 0; beta(isinf(beta)) = 0;
        summation = I*a_n(:,:,end).*alpha - b_n(:,:,end).*beta;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = I*b_n(:,:,end).*alpha - a_n(:,:,end).*beta;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta(:,end);
        
        alpha     = transpose(ric_besselj_derivative(nu,x)).*aux0;
        beta      = transpose(ric_besselj(nu,x)).*aux1;
        alpha(isnan(alpha)) = 0; beta(isinf(beta)) = 0;
        summation = I*a_n(:,:,end).*alpha - b_n(:,:,end).*beta;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = I*b_n(:,:,end).*alpha - a_n(:,:,end).*beta;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta(:,end);        

    elseif (radius(1) <= r) 
        % Outside the multilayer sphere, including the interface
        
        x         = k_m*r;
        
        % Implement (11-239a) in [Balanis1989]    
        alpha     = (transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        alpha(isnan(alpha)) = 0; 
        E_r       =-I*cos(phi)*sum(a_n_scat.*alpha, 2);
        H_r       =-I*sin(phi)*sum(b_n_scat.*alpha, 2)./eta_m;
        
        % Implement (11-239b) in [Balanis1989]    
        alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux1;
        beta      = transpose(ric_besselh(nu,2,x)).*aux0;
        alpha(isnan(alpha)) = 0; beta(isinf(beta)) = 0;
        summation = I*a_n_scat.*alpha - b_n_scat.*beta;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = I*b_n_scat.*alpha - a_n_scat.*beta;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta_m;
        
        % Implement (11-239c) in [Balanis1989]
        alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux0;
        beta      = transpose(ric_besselh(nu,2,x)).*aux1;
        alpha(isnan(alpha)) = 0; beta(isinf(beta)) = 0;
        summation = I*a_n_scat.*alpha - b_n_scat.*beta;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = I*b_n_scat.*alpha - a_n_scat.*beta;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta_m;
    
    else
        % Between the boundary of the inner most sphere and the
        % boundary of the multilayer sphere
        
        idx       = find(r < radius);
        idx       = idx(end);
        
        x         = k_s(:,idx)*r;

        alpha_1   = (transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        alpha_2   = (transpose(ric_besselh_derivative(nu,1,x,2))+transpose(ric_besselh(nu,1,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        alpha_1(isnan(alpha_1)) = 0; alpha_2(isnan(alpha_2)) = 0;
        E_r       = -I*cos(phi)*sum(a_n(:,:,idx).*alpha_1+c_n(:,:,idx).*alpha_2, 2);
        H_r       = -I*sin(phi)*sum(b_n(:,:,idx).*alpha_1+d_n(:,:,idx).*alpha_2, 2)./eta(:,idx);
    
        alpha_1   = transpose(ric_besselh_derivative(nu,2,x)).*aux1;
        alpha_2   = transpose(ric_besselh_derivative(nu,1,x)).*aux1;
        beta_1    = transpose(ric_besselh(nu,2,x)).*aux0;
        beta_2    = transpose(ric_besselh(nu,1,x)).*aux0;
        alpha_1(isnan(alpha_1)) = 0; alpha_2(isnan(alpha_2)) = 0;
        beta_1(isinf(beta_1)) = 0;  beta_2(isinf(beta_2)) = 0;
        summation = I*a_n(:,:,idx).*alpha_1 - b_n(:,:,idx).*beta_1 ...
            +       I*c_n(:,:,idx).*alpha_2 - d_n(:,:,idx).*beta_2;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = I*b_n(:,:,idx).*alpha_1 - a_n(:,:,idx).*beta_1 ...
            +       I*d_n(:,:,idx).*alpha_2 - c_n(:,:,idx).*beta_2;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta(:,idx);
        
        alpha_1   = transpose(ric_besselh_derivative(nu,2,x)).*aux0;
        alpha_2   = transpose(ric_besselh_derivative(nu,1,x)).*aux0;
        beta_1    = transpose(ric_besselh(nu,2,x)).*aux1;
        beta_2    = transpose(ric_besselh(nu,1,x)).*aux1;
        alpha_1(isnan(alpha_1)) = 0; alpha_2(isnan(alpha_2)) = 0;
        beta_1(isinf(beta_1)) = 0;  beta_2(isinf(beta_2)) = 0;
        summation = I*a_n(:,:,idx).*alpha_1 - b_n(:,:,idx).*beta_1 ...
            +       I*c_n(:,:,idx).*alpha_2 - d_n(:,:,idx).*beta_2;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = I*b_n(:,:,idx).*alpha_1 - a_n(:,:,idx).*beta_1 ...
            +       I*d_n(:,:,idx).*alpha_2 - c_n(:,:,idx).*beta_2;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta(:,idx);
    
    end
end


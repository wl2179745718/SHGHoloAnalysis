function  [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                      sphere, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    
    % [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
    %                                              sphere, ...
    %                                              background, ...
    %                                              sensor_location, ...
    %                                              frequency)
    %
    % Calculate the field as a plane wave is scattered by a dielectric
    % sphere centered at the origin. The incident plane wave is
    % polarized in the +x direction propagating in the +z
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
    %   background                   k
    %   --------------------|----------------- 0
    %                       r
    %
    % See Fig. 11-25 in [Balanis1989] for the exact coordinate convention.
    %    
    %
    % Input:
    %
    % radius             Scalar to denote the radius of the sphere (m)
    %
    % sphere             Object of DielectricMaterial
    %
    % background         Object of DielectricMaterial
    %
    % sensor_location    3x1 vector (m)
    %
    % frequency          Nx1 vector (Hz)
    %
    %
    % Output:
    %
    % E_r               Nx1 vector (V/m)
    % E_phi             Nx1 vector (V/m)
    % E_theta           Nx1 vector (V/m)
    % H_r               Nx1 vector (A/m)
    % H_phi             Nx1 vector (A/m)
    % H_theta           Nx1 vector (A/m)

    % See the notes on 2008-05-24 for the coefficients, a_n, b_n, and c_n.

    p = inputParser;    
    p.addRequired('radius', @isnumeric);
    p.addRequired('sphere', @isobject);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParameter('debug', 0, @(x)x==0||x==1);
    p.parse(radius, sphere, background, sensor_location, frequency, varargin{:});
    if (p.Results.debug)
        disp(p.Results);
    end

    EPS_0   = 8.8541878176*1e-12;    
    MU_0    = 4*pi*1e-7;

    nFreq   = length(frequency);

    eta_m   = getIntrinsicImpedance(background, frequency);
    k_m     = getWaveNumber(background, frequency);
    mu_m    = getComplexPermeability(background, frequency)*MU_0;
    eps_m   = getComplexPermittivity(background, frequency)*EPS_0;
    
    eta_s   = getIntrinsicImpedance(sphere, frequency);
    k_s     = getWaveNumber(sphere, frequency);
    mu_s    = getComplexPermeability(sphere, frequency)*MU_0;
    eps_s   = getComplexPermittivity(sphere, frequency)*EPS_0;

    N       = getNMax(p.Results.radius, ...
                      {p.Results.sphere}, ...
                      p.Results.background, ...
                      p.Results.frequency);
    N       = max(N);
    nu      = 1:N;

    [r, theta, phi] = cartToSph(sensor_location(1),sensor_location(2),sensor_location(3));

    % Compute coefficients 
    a_n = 1j.^(-nu).*(2*nu+1)./(nu.*(nu+1)); 
    a_n = ones(nFreq,1)*a_n;

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

    aux1     = ones(nFreq,1)*aux1;
    aux0     = ones(nFreq,1)*aux0;

    iNU = 10;    
    if (p.Results.debug)
        A   = [ric_besselh_derivative(iNU,2,k*radius) -sqrt(eps_m*mu_m)/sqrt(eps_s*mu_s)*ric_besselj_derivative(iNU,k_s*radius);
               ric_besselh(iNU,2,k*radius) -mu_m/mu_s*ric_besselj(iNU,k_s*radius)];
        rhs = -a_n(iNU)*[ric_besselj_derivative(iNU,k*radius); ric_besselj(iNU,k*radius)];
        x   = A\rhs;
        disp(['b_n ' num2str(x(1)) ' d_n ' num2str(x(2))]);

        A   = [ric_besselh(iNU,2,k*radius) -sqrt(eps_m*mu_m)/sqrt(eps_s*mu_s)*ric_besselj(iNU,k_s*radius);
               ric_besselh_derivative(iNU,2,k*radius) -mu_m/mu_s*ric_besselj_derivative(iNU,k_s*radius)];
        rhs = -a_n(iNU)*[ric_besselj(iNU,k*radius); ric_besselj_derivative(iNU,k*radius)];
        x   = A\rhs;
        disp(['c_n ' num2str(x(1)) ' e_n ' num2str(x(2))]);
        disp('------');
    end
    
    
    if (r <= radius)
        num = 1j.*mu_s/sqrt(mu_m)*sqrt(eps_s);
        den =  - sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselh_derivative(nu,2,k_m*radius))...
               + sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselh(nu,2,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius));
        d_n = num*ones(1,N)./den.*a_n;
        
        den = + sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselh(nu,2,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius))...
              - sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselh_derivative(nu,2,k_m*radius));
        e_n = num*ones(1,N)./den.*a_n;
        

        if (p.Results.debug)
            disp(['d_n ' num2str(d_n(iNU)) ' e_n ' num2str(e_n(iNU))]);
            return
        end
        
        x   = k_s*r;
        
        % Implement (11-239a) in [Balanis1989]    
        alpha     = (transpose(ric_besselj_derivative(nu,x,2))+transpose(ric_besselj(nu,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        E_r       = -1j*cos(phi)*sum(d_n.*alpha, 2);
        H_r       = -1j*sin(phi)*sum(e_n.*alpha, 2)./eta_s;
        
        % Implement (11-239b) in [Balanis1989]    
        alpha     = transpose(ric_besselj_derivative(nu,x)).*aux1;
        beta      = transpose(ric_besselj(nu,x)).*aux0;
        summation = 1j*d_n.*alpha - e_n.*beta;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = 1j*e_n.*alpha - d_n.*beta;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta_s;
        
        % Implement (11-239c) in [Balanis1989]
        alpha     = transpose(ric_besselj_derivative(nu,x)).*aux0;
        beta      = transpose(ric_besselj(nu,x)).*aux1;
        summation = 1j*d_n.*alpha - e_n.*beta;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = 1j*e_n.*alpha - d_n.*beta;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta_s;        
    else
        
        num =  + sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselj(nu,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius)) ...
               - sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselj_derivative(nu,k_m*radius));
        den =  + sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselh_derivative(nu,2,k_m*radius))...
               - sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselh(nu,2,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius));
        b_n = num./den.*a_n;
        
        num = + sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselj_derivative(nu,k_m*radius))...
              - sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius));
        den = + sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselh(nu,2,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius))...
              - sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselh_derivative(nu,2,k_m*radius));
        c_n = num./den.*a_n;
        
        if (p.Results.debug)
            disp(['b_n ' num2str(b_n(iNU)) ' c_n ' num2str(c_n(iNU))]);
            return
        end
        
        x   = k_m*r;
        
        % Implement (11-239a) in [Balanis1989]    
        alpha     = (transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
            .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
        E_r       = -1j*cos(phi)*sum(b_n.*alpha, 2);
        H_r       = -1j*sin(phi)*sum(c_n.*alpha, 2)./eta_m;
        
        % Implement (11-239b) in [Balanis1989]    
        alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux1;
        beta      = transpose(ric_besselh(nu,2,x)).*aux0;
        summation = 1j*b_n.*alpha - c_n.*beta;
        E_theta   = cos(phi)./x.*sum(summation,2);
        summation = 1j*c_n.*alpha - b_n.*beta;
        H_theta   = sin(phi)./x.*sum(summation,2)./eta_m;
        
        % Implement (11-239c) in [Balanis1989]
        alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux0;
        beta      = transpose(ric_besselh(nu,2,x)).*aux1;
        summation = 1j*b_n.*alpha - c_n.*beta;
        E_phi     = sin(phi)./x.*sum(summation,2);
        summation = 1j*c_n.*alpha - b_n.*beta;
        H_phi     =-cos(phi)./x.*sum(summation,2)./eta_m;
    end
    
end


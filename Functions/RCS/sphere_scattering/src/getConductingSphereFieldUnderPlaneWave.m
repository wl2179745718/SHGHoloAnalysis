function [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getConductingSphereFieldUnderPlaneWave(radius, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
    
    % [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = getConductingSphereFieldUnderPlaneWave(radius, ...
    %                                        background, ...
    %                                        sensor_location, ...
    %                                        background, ...
    %                                        frequency)
    % Calculate the field as a plane wave is scattered by a PEC sphere
    % centered at the origin.  The incident plane wave is polarized in
    % the +x direction propagating in the +z direction.
    %
    % See Fig. 11-25 in [Balanis1989] for the exact coordinate convention.
    %    
    %
    % Input:
    %
    % radius             Scalar to denote the radius of the sphere (m)
    %
    % background         Object of DielectricMaterial
    %
    % sensor_location    3x1 vector in (m)
    %
    % frequency          Nx1 vector in (Hz)
    % 
    % Output:
    %
    % E_r               Nx1 vector (V/m)
    % E_phi             Nx1 vector (V/m)
    % E_theta           Nx1 vector (V/m)
    % H_r               Nx1 vector (A/m)
    % H_phi             Nx1 vector (A/m)
    % H_theta           Nx1 vector (A/m)
    
    % The coefficients, a_n, b_n, and c_n are in (11-231a), (11-238a), and
    % (11-238b) in [Balanis1989].  And the scattered field is in
    % (11-239).

    p = inputParser;    
    p.addRequired('radius', @isnumeric);
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParameter('debug', 0, @(x)x==0||x==1);
    p.parse(radius, background, sensor_location, frequency, varargin{:});
    if (p.Results.debug)
        disp(p.Results);
    end

    nFreq   = length(frequency);
    omega   = 2*pi*reshape(frequency,[nFreq 1]);
    eta     = getIntrinsicImpedance(background, frequency);
    k       = getWaveNumber(background, frequency);

    N       = 50; % Order of the Bessel(Hankel). Important. If this number is
                  % not large enough, the series does not converge.
    nu      = 1:N;

    [r, theta, phi] = cartToSph(sensor_location(1),sensor_location(2),sensor_location(3));

    if (r < radius)
        % Field inside the sphere is 0
        E_r = 0;
        E_theta = 0;
        E_phi = 0;
        H_r = 0;
        H_theta = 0;
        H_phi = 0;
        return;
    end

    % Compute coefficients as in (11-231a), (11-238a), (11-239b)
    a_n = 1j.^(-nu).*(2*nu+1)./(nu.*(nu+1)); 
    a_n = ones(nFreq,1)*a_n;
    b_n = -a_n.*transpose(ric_besselj_derivative(nu,k*radius))./transpose(ric_besselh_derivative(nu,2,k*radius));
    c_n = -a_n.*transpose(ric_besselj(nu,k*radius))./transpose(ric_besselh(nu,2, k*radius));
    
    % aux0 denote the expression assoc_legendre(nu,1,cos(theta))/sin(theta). Here I
    % am using an recursive relation to compute aux0, which avoids the
    % numerical difficulty when theta == 0 or PI.
    aux0    = zeros(1, length(nu)); 
    aux0(1) = -1; 
    aux0(2) = -3*cos(theta);  
    for n = 2:nu(end-1) 
        aux0(n+1) = (2*n+1)/n*cos(theta)*aux0(n) - (n+1)/n*aux0(n-1);
    end
    
    % aux1 denote the expression
    % sin(theta)*assoc_legendre_derivative(nu,1,cos(theta)).  Here I am
    % also using an recursive relation to compute aux1 from aux0,
    % which avoids numerical difficulty when theta == 0 or PI.
    aux1    = zeros(1, length(nu));
    aux1(1) = cos(theta); 
    for n = 2:nu(end)
        aux1(n) = (n+1)*aux0(n-1)-n*cos(theta)*aux0(n);
    end

    aux1     = ones(nFreq,1)*aux1;
    aux0     = ones(nFreq,1)*aux0;
    
    x   = k*r;
    
    % Implement (11-239a) in [Balanis1989]    
    alpha     = ( transpose(ric_besselh_derivative(nu,2,x,2))+transpose(ric_besselh(nu,2,x)))...
        .*transpose(assoc_legendre(nu,1,cos(theta))*ones(1,nFreq));
    E_r       = -1j*cos(phi)*sum(b_n.*alpha, 2);
    H_r       = -1j*sin(phi)*sum(c_n.*alpha, 2)./eta;
    
    % Implement (11-239b) in [Balanis1989]    
    alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux1;
    beta      = transpose(ric_besselh(nu,2,x)).*aux0;
    summation = 1j*b_n.*alpha - c_n.*beta;
    E_theta   = cos(phi)./x.*sum(summation,2);
    summation = 1j*c_n.*alpha - b_n.*beta;
    H_theta   = sin(phi)./x.*sum(summation,2)./eta;
    
    % Implement (11-239c) in [Balanis1989]
    alpha     = transpose(ric_besselh_derivative(nu,2,x)).*aux0;
    beta      = transpose(ric_besselh(nu,2,x)).*aux1;
    summation = 1j*b_n.*alpha - c_n.*beta;
    E_phi     = sin(phi)./x.*sum(summation,2);
    summation = 1j*c_n.*alpha - b_n.*beta;
    H_phi     =-cos(phi)./x.*sum(summation,2)./eta;

end

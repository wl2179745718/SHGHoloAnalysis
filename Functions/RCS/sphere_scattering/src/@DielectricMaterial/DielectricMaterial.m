classdef DielectricMaterial
    properties (SetAccess = private)
        epsilon_infty = 1;
        sigma_e = 0.0;
        mu_r = 1;
        sigma_m = 0.0;
        th_conv = 'engr';
    end
    
    methods
        function this = DielectricMaterial(varargin)
        % this = DielectricMaterial(epsilon_infty, sigma_e, mu_r, sigma_m)
        %
        % this = DielectricMaterial(epsilon_infty,sigma_e) assumes a non-magnetic
        % material: mu_r =1 and sigma_m = 0
        %
            if (nargin == 2)
                this.epsilon_infty = cell2mat(varargin(1));
                this.sigma_e       = cell2mat(varargin(2));
                this.mu_r          = 1;
                this.sigma_m       = 0;
                this.th_conv       = 'engr';
            elseif (nargin ==4)
                this.epsilon_infty = cell2mat(varargin(1));
                this.sigma_e       = cell2mat(varargin(2));
                this.mu_r          = cell2mat(varargin(3));
                this.sigma_m       = cell2mat(varargin(4));
                this.th_conv       = 'engr';
            elseif (nargin ==5)
                this.epsilon_infty = cell2mat(varargin(1));
                this.sigma_e       = cell2mat(varargin(2));
                this.mu_r          = cell2mat(varargin(3));
                this.sigma_m       = cell2mat(varargin(4));
                this.th_conv       = cell2mat(varargin(5));
            else
                error('??? Invalid argument to constructor.');
            end
        end%DielectricMaterial
    end%methods
    
    methods
        function z = convertToAbsorptionDepth(this, distance, frequency)
        % z = convertToAbsorptionDepth(this, distance, frequency)
        %
            
            z = getAbsorptionDepth(this, frequency);
            z = distance./z;
        end
        
        function z = convertToLength(this, f, x)
        % z = convertToLength(this, f, x) 
        %
        % f    frequency
        % x    length (skin depth)
            
            z = x.*getAbsorptionDepth(this, f);
        end
        
        function x = getAbsorptionDepth(this, frequency)
        % x = getAbsorptionDepth(this, frequency) calculates the absorption
        % depth from the wave number. The absorption depth is always
        % positive regardless the choice of the time-harmonic factor.
            
            k = getWaveNumber(this, frequency);
            x = abs(1./imag(k));
        end

        function alpha = getAttenuationConstant(this, frequency, varargin)
        % alpha = getAttenuationConstant(this, frequency, varargin)
        %
        % alpha = getAttenuationConstant(this, frequency, 'dB') gives the
        % value in 'dB/m'
        %

            p = inputParser;
            p.addOptional('scale', 1, @(x)strcmp(x, 'dB'))
            p.parse(varargin{:});
            
            k = getWaveNumber(this, frequency);
            if (p.Results.scale == 'dB')
                alpha = -20*log10(exp(1)).*imag(k);
            else
                alpha = imag(k);
            end
            
        end            

        function mu_r = getComplexPermeability(this, frequency)
        % mu_r = getComplexPermeability(this, frequency) computes the relative
        % complex permeability.
        %
        % Input:
        %
        % frequency    Nx1 vector (Hz)
        %
            
            MU_0      = 4*pi*1e-7;
            mu_r      = this.mu_r + this.sigma_m./(1j*2*pi*frequency*MU_0);
            mu_r(frequency==0) = this.mu_r;
            
            
            
        end

        function epsilon_r = getComplexPermittivity(this, frequency)
        % epsilon_r = getComplexPermittivity(this, frequency) computes the
        % relative complex permittivity.
        %
        % Input:
        %
        %  frequency   Nx1 vector (Hz)
        %
        %
            
            EPS_0          = 8.8541878176e-12;    
            idx            = find(frequency==0);

            switch this.th_conv
              case 'engr'
                epsilon_r      = this.epsilon_infty + this.sigma_e./(1j*2*pi*frequency*EPS_0);
                epsilon_r(idx) = this.epsilon_infty;
              case 'phys'
                epsilon_r      = this.epsilon_infty + this.sigma_e./((-1j)*2*pi*frequency*EPS_0);
                epsilon_r(idx) = this.epsilon_infty;
              otherwise
                error('Invalid');
            end

        end

        function ref_idx = getComplexRefractiveIndex(this, frequency)
        % ref_idx = getComplexRefractiveIndex(this, frequency) computes the
        % refractive index.
        %
        % Input:
        %
        % frequency    Nx1 vector (Hz)
        %
            
            eps_r = getComplexPermittivity(this, frequency);
            mu_r  = getComplexPermeability(this, frequency);
            ref_idx =  sqrt(eps_r.*mu_r);
            

        end

        function [v_g, f] = getGroupVelocity(this, frequency)
        % v_g = getGroupVelocity(this, frequency) evalutes the group velocity
        % by numerically differentiating the angular frequency with respect to
        % the wave number.
        %
        %

            omega = 2*pi*frequency;
            k     = getWaveNumber(this, frequency);
            v_g   = diff(omega)./diff(real(k));
            f     = (frequency(1:end-1)+frequency(2:end))/2;

        end

        function eta = getIntrinsicImpedance(this, frequency)
        % eta = getIntrinsicImpedance(this, frequency)
        %
        % Input:
        % 
        % frequency     Nx1 vector (Hz)
        %    

            EPS_0 = 8.8541878176e-12;
            MU_0  = 4*pi*1e-7;
            ETA_0 = sqrt(MU_0/EPS_0);
            eta   = ETA_0*sqrt(getComplexPermeability(this, frequency)./ ...
                               getComplexPermittivity(this, frequency));
            
            
        end

        function [mu_r, sigma_m] = getPermeability(this)
        % [mu_r, sigma_m] = getPermeability(this, frequency)
        %
        % Input:
        % 
        % frequency    Nx1 vector (Hz)
            
            
            mu_r    = this.mu_r;
            sigma_m = this.sigma_m;
        end

        function [epsilon_r, sigma_e] = getPermittivity(this, frequency)
        % [epsilon_r, sigma_e] = getPermittivity(this)
        %
        % [epsilon_r, sigma_e] = getPermittivity(this, frequency)
        %
        % Input:
        % 
        % frequency     Nx1 vector (Hz)
        %
            
            epsilon_r = ones(size(frequency))*this.epsilon_infty;
            sigma_e   = ones(size(frequency))*this.sigma_e;
            
            
            
        end

        function v_p = getPhaseVelocity(this, frequency)
        % v_p = getPhaseVelocity(this, frequency)

            omega = 2*pi*frequency;
            k     = getWaveNumber(this, frequency);
            v_p   = omega./real(k);
            
            
        end

        function epsilon_r = getSymbolicPermittivity(this)
        %

            EPSILON_0 = 8.8541878176e-12;
            
            I = sqrt(-1);
            syms omega
            
            switch this.th_conv
              case 'engr'
                epsilon_r = this.epsilon_infty + this.sigma_e/(+I*omega*EPSILON_0);
              case 'phys'
                epsilon_r = this.epsilon_infty + this.sigma_e/(-I*omega*EPSILON_0);
              otherwise
                error('Error');
            end
            
            
        end

        function lambda = getWavelength(this, frequency)
        % lambda = getWaveLength(this, frequency)

            k      = getWaveNumber(this, frequency);
            lambda = 2*pi./real(k);
            
            
        end

        function k = getWaveNumber(this, frequency)
        % k = getWaveNumber(this, frequency)
        %
        % frequency can be a Nx1 vector
        % wave_number may be complex
        %

            C_0          = 2.997924580003452e+08;
            permittivity = getComplexPermittivity(this, frequency);
            permeability = getComplexPermeability(this, frequency);
            k            = 2*pi*frequency.*sqrt(permittivity.*permeability)/C_0;
        end
        
    end%methods
end%classdef

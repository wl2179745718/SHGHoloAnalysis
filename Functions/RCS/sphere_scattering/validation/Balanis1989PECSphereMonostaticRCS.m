function Balanis1989PECSphereMonostaticRCS(varargin)
% Duplicate Fig. 11-26 in [Balanis1989].
%
%
    addpath('../src');
    
    radius          = 1.0; %m
    ratio           = (0.01:0.01:1.6)';
    wave_length     = radius./ratio;
    background      = DielectricMaterial(1, 0);
    frequency       = getPhaseVelocity(background, 1)./wave_length;
    sensor_location = [0; 0; -max(wave_length)*20];

    % Set epsilon_r is an extremely large value and mu_r an extremely
    % small value, while maintain their product to 1. Such a sphere
    % has the same boundry condition as a PEC sphere. See
    % [Harrington2001].
    sphere          = DielectricMaterial(1e8, 0, 1e-8, 0); 
    [E_r, E_theta, E_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                      sphere, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency);
    E               = [E_r, E_theta, E_phi];
    mono_RCS        = 4*pi*norm(sensor_location)^2*sum(E.*conj(E),2);
    y               = mono_RCS/(pi*radius^2);
    
    figure;
    semilogy(ratio, y, 'LineWidth', 2); grid('on');
    set(gca,'fontsize', 12); set(gca,'xlim',[0 1.6], 'ylim', [0.01, 5]);
    set(gca,'ytick',[0.01 0.05 0.1 0.5 1 5])
    xlabel('Sphere radius in wavelengths');
    ylabel('Normalized monostatic RCS');
end


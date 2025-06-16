function Djordjevic2004()
% Plot the monostatic RCS of a conducting sphere as in Fig. 8 and
% plot the monostatic RCS of a dielectric sphere as in Fig. 9 in
% [Djordjevic2004]. The relative permittivity of the dielectric
% sphere is 2.25.

    addpath('../src');

    radius          = 1.0; % m
    ratio           = (0.01:0.006:2.6)';
    wavelength      = radius./ratio;
    background      = DielectricMaterial(1, 0);
    frequency       = getPhaseVelocity(background, 1)./wavelength;
    sensor_location = [0; 0; -max(wavelength)*20]; % 20 wavelength away

    % The conducting sphere can be effectively mimicked by a material with
    % a very large permittivity and a very small permeability.
    sphere          = DielectricMaterial(1e8, 0, 1e-8, 0);

    [E_r,E_theta,E_phi] = getMultilayerSphereFieldUnderPlaneWave(radius, ...
                                                      {sphere}, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency);
    E               = [E_r E_theta E_phi];
    mono_RCS        = 4*pi*norm(sensor_location)^2*sum(E.*conj(E),2);
    y               = mono_RCS/(pi*radius^2);

    set(figure, 'color', 'white');
    semilogy(ratio, y, 'LineWidth', 2); 
    grid('on');
    set(gca,'fontsize', 12); 
    set(gca,'xlim',[0 2.0], 'ylim', [0.01, 100]);
    set(gca,'xtick', 0:0.2:2.0, 'ytick',[0.01 0.1 1 10 100])
    xlabel('Sphere radius / free-space wavelengths');
    ylabel('Normalized RCS of a conducting sphere');


    % The dielectric sphere
    sphere              = DielectricMaterial(2.25,0);

    % using getDielectricSphereFieldUnderPlaneWave()
    [E_r,E_theta,E_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                      sphere, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency);
    E        = [E_r E_theta E_phi];
    mono_RCS = 4*pi*norm(sensor_location)^2*sum(E.*conj(E),2);
    y1       = mono_RCS/(pi*radius^2);

    % using getMultilayerSphereFieldUnderPlaneWave()
    [E_r,E_theta,E_phi] = getMultilayerSphereFieldUnderPlaneWave(radius, ...
                                                      {sphere}, ...
                                                      background, ...
                                                      sensor_location, ...
                                                      frequency);
    E        = [E_r E_theta E_phi];
    mono_RCS = 4*pi*norm(sensor_location)^2*sum(E.*conj(E),2);
    y2       = mono_RCS/(pi*radius^2);

    set(figure, 'color', 'white');
    semilogy(ratio, y1, 'LineWidth', 2); 
    hold('on');
    semilogy(ratio, y2, 'r--', 'LineWidth', 2); 
    grid('on');
    set(gca,'fontsize', 12); 
    set(gca,'xlim',[0 2.0], 'ylim', [0.01, 100]);
    set(gca,'xtick', 0:0.2:2.0, 'ytick',[0.01 0.1 1 10 100])
    xlabel('Sphere radius / free-space wavelength');
    ylabel('Normalized RCS of a dielectric sphere');
    legend('getDielectricSphere()', 'getMultilayerSphere()');



function Ruck1970PECSphereBistaticRCS(varargin)
% See Fig. 3-1 on pp. 141 in [Ruck1970] about the geometry. Notice that
% he assumes a plane wave propagating in the -z direction, which is
% opposite to my assumption.
%
% Plot the bistatic RCS of a conducting sphere, which is shown in
% Fig. 3-5 on pp. 152 in [Ruck1970]. This function calls one of the
% conducting, dielectric, and multilayer sphere code. They give
% identical results. See notes on 2010-11-08.

    addpath('../src');    
    
    p = inputParser;
    p.addParameter('code', 'conducting', @(x)strcmp(x,'conducting') || strcmp(x,'die') || strcmp(x, 'multi'));
    p.addParameter('debug', 0, @(x)x==0||x==1);
    p.parse(varargin{:});
    if (p.Results.debug)
        disp(p.Results.debug);
    end
    
    k_oA            = 7.7;
    background      = DielectricMaterial(1, 0);
    frequency       = 1e9;
    k_o             = getWaveNumber(background, frequency);
    radius          = k_oA/k_o;
    sphere          = DielectricMaterial(1e5, 0, 1e-5, 0);
    theta           = linspace(pi, 0, 2^6)';
    
    sensor_distance = 10000*2*pi/k_o;
    E               = zeros(length(theta), 3);
    H               = zeros(length(theta), 3);
    for i = 1:length(theta)
        sensor_location = [0; 
                           sensor_distance*sin(theta(i)); 
                           sensor_distance*cos(theta(i))];
        switch p.Results.code
          case 'conducting'
            [E(i,1),E(i,2),E(i,3),H(i,1),H(i,2),H(i,3)] = getConductingSphereFieldUnderPlaneWave(radius, ...
                                                              background, ...
                                                              sensor_location, ...
                                                              frequency);
          case 'die'
            [E(i,1),E(i,2),E(i,3),H(i,1),H(i,2),H(i,3)] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                              sphere, ...
                                                              background, ...
                                                              sensor_location, ...
                                                              frequency);
          case 'multi'
            [E(i,1),E(i,2),E(i,3),H(i,1),H(i,2),H(i,3)] = getMultilayerSphereFieldUnderPlaneWave(radius, ...
                                                              {sphere}, ...
                                                              background, ...
                                                              sensor_location, ...
                                                              frequency);
        end
    end%for
    
    % RCS definition in Fig. 11-26 in [Balanis1989]
    RCS = 4*pi*sensor_distance^2*sum(E.*conj(E),2);
    RCS = RCS./(pi*radius.^2);

    figure; 
    semilogy((180/pi)*theta, RCS,'b-','linewidth', 2);
    hold('on');
    
    grid('on');
    set(gca, 'xlim', [0 180]);
    set(gca, 'xtick', 0:30:180,'plotboxaspectratio',[1 0.5 1]);
    set(gca, 'fontsize', 12, 'fontname', 'arial');
    xlabel('\theta (deg)');
    ylabel('Normalized bi-static RCS');

end


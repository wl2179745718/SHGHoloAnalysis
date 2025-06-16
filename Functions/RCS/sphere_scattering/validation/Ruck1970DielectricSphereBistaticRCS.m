function Ruck1970DielectricSphereBistaticRCS(varargin)
% plot the radar cross-section of a dielectric sphere with a
% refractive index of 1.61, which is fig. 3-7 on pp. 162 of
% [ruck1970].
%

    addpath('../src');
    
    k_oa       = 3.0;
    background = DielectricMaterial(1, 0);
    frequency  = 1e9;
    k_o        = getWaveNumber(background, frequency);
    radius     = k_oa/k_o;
    sphere     = DielectricMaterial(1.61^2, 0.0);
    theta      = linspace(pi, 0, 2^6)';
    
    sensor_distance = 3*2*pi/k_o; 

    for i = 1:length(theta)
        sensor_location = [0;
                           sensor_distance*sin(theta(i));
                           sensor_distance*cos(theta(i))];
        
        [e_r, e_theta, e_phi] = getDielectricSphereFieldUnderPlaneWave(radius, ...
                                                          sphere, ...
                                                          background, ...
                                                          sensor_location, ...
                                                          frequency);
        
        e = [e_r; e_theta; e_phi];
        % disp(e);
        
        if (theta(i)==0) 
            rcs_e(i,1) = 4*pi*sensor_distance^2*(e(2)'*e(2)); 
        else
            rcs_e(i,1) = 4*pi*sensor_distance^2*(e(3)'*e(3)); 
        end
    end
    
    % normalize the rcs
    rcs_e = rcs_e./(pi*radius.^2);
    
    figure;
    set(gcf, 'color', 'white');
    semilogy(theta/pi*180, rcs_e, 'b-', 'linewidth',2); 
    hold('on');
    grid('on');
    set(gca, 'xlim', [0, 180]);
    set(gca,'xtick', 0:30:180);
    xlabel('\theta (deg)'); 
    ylabel('Normalized bistatic rcs - e');
end



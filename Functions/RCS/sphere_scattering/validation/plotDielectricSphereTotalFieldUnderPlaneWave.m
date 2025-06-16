function plotDielectricSphereTotalFieldUnderPlaneWave(varargin)
% Plot the total field along the three transects, as a plane wave is
% scattered by a sphere centered at the origin. The plane wave
% propagates in the +z direction. See Fig. 11-25 in [Balanis1989] for
% the geometry.
%
%
% Optional:
%
% 'radius'        radius of the sphere relative to the background
%                 wavelength). Default to 1.
%
% 'sphere'        object of DielectricMaterial. Default to
%                 DielectricMaterial(4.0, 0.0).
%
% 'background'    object of DielectricMaterial. Default to
%                 DielectricMaterial(1.0, 0.0).
%
% 'frequency'     scalar to denote the frequency of the plane
%                 wave. Default to 1 Hz.
%
% 'domain_size'   6x1 vector to specify the size of the domain in
%                 the format  [x_min; x_max; y_min; y_max; z_min; z_max]
%                 (relative to the background wavelength). Default
%                 to [-0.2; 2; -0.2; 2; -0.2; 2].
%
% 'N'             3x1 vector to specify the number points to
%                 evaluate the fields along the three axis.
%                 Default to [44; 44; 44].
%
% 'S'             3x1 vector to specify the slice on which the
%                 fields are evaluated. Default to [4;4;60]
%
% 'workspace_name'  workspace name to store the
%                   data. Default to 'dielectric_sphere'.
%
    
    addpath('../src');
    
    p = inputParser;
    p.addParameter('radius', 1.0, @isnumeric);
    p.addParameter('sphere', DielectricMaterial(4.0, 0.0), @isobject);
    p.addParameter('background', DielectricMaterial(1.0, 0.0), @isobject);
    p.addParameter('frequency', 1.0, @isnumeric);
    p.addParameter('domain_size', [-0.2; +2; -0.2; +2; -0.2; +2], @isvector);
    p.addParameter('N', [45; 45; 45], @isvector);
    p.addParameter('S', [5;   5;  5], @isvector);
    p.addParameter('workspace_name', 'dielectric_sphere', @ischar);
    p.addParameter('verbose', 0, @(x)x==0||x==1);
    p.parse(varargin{:});
    if (p.Results.verbose)
        disp(p.Results);
    end

    wavelength = getWavelength(p.Results.background, ...
                               p.Results.frequency);
    x          = linspace(p.Results.domain_size(1), p.Results.domain_size(2), p.Results.N(1));
    y          = linspace(p.Results.domain_size(3), p.Results.domain_size(4), p.Results.N(2));
    z          = linspace(p.Results.domain_size(5), p.Results.domain_size(6), p.Results.N(3));
    
    E_r        = zeros([length(x) length(y) length(z)]);
    E_theta    = zeros([length(x) length(y) length(z)]);
    E_phi      = zeros([length(x) length(y) length(z)]);
    H_r        = zeros([length(x) length(y) length(z)]);
    H_theta    = zeros([length(x) length(y) length(z)]);
    H_phi      = zeros([length(x) length(y) length(z)]);

    for iX = 1:p.Results.N(1)
        for iY = 1:p.Results.N(2)
            for iZ = 1:p.Results.N(3)
                if (iX==p.Results.S(1) ...
                    || iY==p.Results.S(2) ...
                    || iZ==p.Results.S(3))
                    [E_r(iX,iY,iZ), E_theta(iX,iY,iZ), E_phi(iX,iY,iZ), ...
                     H_r(iX,iY,iZ), H_theta(iX,iY,iZ), H_phi(iX,iY,iZ)] ...
                        = getDielectricSphereFieldUnderPlaneWave(p.Results.radius*wavelength, ...
                                                                 p.Results.sphere, ...
                                                                 p.Results.background, ...
                                                                 [x(iX);y(iY);z(iZ)]*wavelength,...
                                                                 p.Results.frequency);
                end                
            end
        end
    end

    for iX = 1:p.Results.N(1)
        for iY = 1:p.Results.N(2)
            for iZ = 1:p.Results.N(3)
                if (iX==p.Results.S(1) ...
                    || iY==p.Results.S(2) ...
                    || iZ==p.Results.S(3))
                    sensor_location = [x(iX); y(iY); z(iZ)];
                    if (norm(sensor_location) > p.Results.radius) % outside the sphere
                        [E_x, H_y]      = getPlaneWaveUsingCartesianExpansion(p.Results.background, ...
                                                                          sensor_location*wavelength, ...
                                                                          p.Results.frequency);
                        
                        % Vector transformation from the cartesian to spherical
                        % coordinate. See (II-12a) on pg. 924 in
                        % [Balanis1989].
                        [~, theta, phi] = cartToSph(x(iX), y(iY), z(iZ));
                        A               = [+sin(theta)*cos(phi) sin(theta)*sin(phi) +cos(theta);
                                           +cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta);
                                           -sin(phi)            cos(phi)             0];
                        temp              = A*[E_x;0;0];
                        E_r(iX,iY,iZ)     = E_r(iX,iY,iZ)    +temp(1);
                        E_theta(iX,iY,iZ) = E_theta(iX,iY,iZ)+temp(2);
                        E_phi(iX,iY,iZ)   = E_phi(iX,iY,iZ)  +temp(3);
                        temp              = A*[0;H_y;0];
                        H_r(iX,iY,iZ)     = H_r(iX,iY,iZ)    +temp(1);
                        H_theta(iX,iY,iZ) = H_theta(iX,iY,iZ)+temp(2);
                        H_phi(iX,iY,iZ)   = H_phi(iX,iY,iZ)  +temp(3);
                    end
                end
            end
        end
    end
    
    if (~isempty(p.Results.workspace_name))
        save(p.Results.workspace_name, 'p', 'x', 'y', 'z', ...
             'E_r', 'E_theta', 'E_phi', ...
             'H_r', 'H_theta', 'H_phi');
    end
    
    set(figure,'color','white'); 
    plotSphere(x,y,z,E_r,p.Results.radius,p.Results.S);
    title('E_r');
    set(figure,'color','white'); plotSphere(x,y,z,E_theta,p.Results.radius,p.Results.S);
    title('E_\theta');
    set(figure,'color','white'); plotSphere(x,y,z,E_phi,p.Results.radius,p.Results.S);
    title('E_\phi');
    set(figure,'color','white'); plotSphere(x,y,z,H_r,p.Results.radius,p.Results.S);
    title('H_r');
    set(figure,'color','white'); plotSphere(x,y,z,H_theta,p.Results.radius,p.Results.S);
    title('H_\theta');
    set(figure,'color','white'); plotSphere(x,y,z,H_phi,p.Results.radius,p.Results.S);
    title('H_\phi');
end

function plotSphere(x,y,z,data,radius,slice_idx)
% x      Rx1 vector
% y      Sx1 vector
% z      Tx1 vector
% data   RxSxT matrix
%

    temp_x = [min(x);   max(x)];
    temp_y = [0; 0];
    temp_z = [0; 0];
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    hold('on');
    temp_x = [0; 0];
    temp_y = [min(y);   max(y)];
    temp_z = [0; 0];
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    temp_x = [0; 0];
    temp_y = [0; 0];
    temp_z = [min(z);   max(z)];
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    
    temp_x = zeros(64,1);
    temp_y = radius*sin(linspace(0,pi/2,64));
    temp_z = radius*cos(linspace(0,pi/2,64));
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    
    temp_x = radius*sin(linspace(0,pi/2,64));
    temp_y = zeros(64,1);
    temp_z = radius*cos(linspace(0,pi/2,64));
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    
    temp_x = radius*sin(linspace(0,pi/2,64));
    temp_y = radius*cos(linspace(0,pi/2,64));
    temp_z = zeros(64,1);
    plot3(temp_x,temp_y,temp_z,'k-','linewidth',1.5);
    
    
    % shuffle the data
    [yi, xi, zi] = meshgrid(x, y, z);
    data      = permute(data,[2 1 3]);
    h         = slice(yi,xi,zi,abs(data), ...
                      x(slice_idx(1)), y(slice_idx(2)), z(slice_idx(3)),'nearest');
    
    axis('equal', 'tight'); 
    set(h,'linestyle','none','facealpha',0.7)
    % temp = abs(data(:));
    % caxis([min(temp) max(temp)]); 
    colorbar('EastOutside'); 

    % set the right view
    view([135 45]);
    
    set(gca,'fontname', 'arial', 'fontsize',16)
    xlabel('x (\lambda_0)');
    ylabel('y (\lambda_0)');
    zlabel('z (\lambda_0)');
end

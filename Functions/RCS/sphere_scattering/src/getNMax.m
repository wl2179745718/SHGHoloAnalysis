function N_max = getNMax(radius, sphere, background, frequency)
% N_max = getNMax(radius, sphere, background, frequency) implements
% [Wiscombe1980] to determine the number of terms in the series.
%
% Input:
%
% radius       N x 1 vector in the descending order
% sphere       N x 1 DielectricMaterial/DebyeMaterial cell array in the
%              same order as radius
% background   DielectricMaterial
% frequency    M x 1 vector
%
% Output:
%
% N_max        M x 1 vector
%
    
    radius    = reshape(radius, 1, length(radius));
    sphere    = reshape(sphere, 1, length(sphere));
    frequency = reshape(frequency, length(frequency), 1);
    
    % Follow the index convention in [Yang2003]. The sphere is
    % specified from the inner most to the outter most.
    radius    = fliplr(radius);
    sphere    = fliplr(sphere);
    
    k_m       = getWaveNumber(background, frequency);
    x         = k_m*radius; % x is a matrix
    x         = abs(x);     % x may be complex the background is lossy
    
    N_m       = getComplexRefractiveIndex(background, frequency);
    m         = zeros(length(frequency), length(sphere));
    for iSphere = 1:length(sphere)
        % relative refractive index
        m(:, iSphere)  = getComplexRefractiveIndex(sphere{iSphere}, frequency)./N_m;
    end
    
    N_max      = ones(length(frequency), 1);
    for iFreq = 1:length(frequency)
        % if (0.02 <= x(iFreq, end) & x(iFreq, end) < 8)
        if (x(iFreq, end) < 8) 
            % this is to be consistent with scattnlay. See Line 25 in nmie.c
            N_stop = x(iFreq, end)+4   *x(iFreq, end).^(1/3)+1; 
        elseif (8.00 <= x(iFreq, end) && x(iFreq, end) < 4200)
            N_stop = x(iFreq, end)+4.05*x(iFreq, end).^(1/3)+2;
        else 
            % 4200 <= x(iFreq, end) & x(iFreq, end) < 20000);
            N_stop = x(iFreq, end)+4   *x(iFreq, end).^(1/3)+2;
        end
        
        N_max(iFreq) = max([N_stop ...
                            abs(m(iFreq,:).*x(iFreq,:)) ...
                            abs(m(iFreq,2:end).*x(iFreq,1:end-1))])+15;
        N_max(iFreq) = round(N_max(iFreq));
    end
    
    
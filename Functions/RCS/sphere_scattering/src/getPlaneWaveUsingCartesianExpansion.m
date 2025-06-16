function [Ex, Hy] = getPlaneWaveUsingCartesianExpansion(background, ...
                                                      sensor_location, ...
                                                      frequency, ...
                                                      varargin)
% [Ex, Hy] = getPlaneWaveUsingCartesianExpansion(background,
%                                                sensor_location, 
%                                                frequency) 
%
% Calculate a planewave assuming it propagating in the +z
% direction polarized in the +x direction. Uses e^{j\omega t} as
% the time-harmonic factor.
% 
% Input:
%
% background          DielectricMaterial
% sensor_location     [x;y;z] (m)
% frequency           Nx1 vector (Hz)
%
% Optional:
%
% 'time'              scalar (s). Computes the field at a
%                     particular time instance. Default to [].
%
% Output:
%
% Ex                  Nx1 vector (V/m)
% Hy                  Nx1 vector (A/m)


    p = inputParser;
    p.addRequired('background', @isobject);
    p.addRequired('sensor_location', @isvector);
    p.addRequired('frequency', @isnumeric);
    p.addParameter('time', [], @isnumeric);
    p.parse(background, sensor_location, frequency, varargin{:});
    
    frequency   = reshape(frequency, length(frequency), 1);
    wave_number = getWaveNumber(background, frequency);
    Ex          = exp(-1j*wave_number*sensor_location(3));
    Hy          = +1./getIntrinsicImpedance(background, frequency).*Ex;

    
    if (~isempty(p.Results.time))
        omega = 2*pi*frequency;
        Ex    = real(Ex*exp(1j*omega*p.Results.time));
        Hy    = real(Hy*exp(1j*omega*p.Results.time));
    end
    
    


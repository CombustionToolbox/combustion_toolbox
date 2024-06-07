classdef Constants < handle
    % Class with constants data
    %
    % Attributes:
    %     release (char): Release of the Combustion Toolbox
    %     date (char): Date of the release
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     G (float): Standard gravity [m/s2]
    %
    % Examples:
    %     * R0 = Constants.R0
    %     * g = Constants.G
    %     * release = Constants.release

    
    properties (Constant)
        R0      = 8.31446261815324 % Universal gas constant [J/(K mol)]
        G       = 9.80665          % Standard gravity [m/s2]
        NA      = 6.0221415e23     % Avogadro's number [molecule/mol]
        release = 'v1.1.0'         % Release version
        date    = '07 Jun 2024'    % Release date
    end
    
end
classdef Constants < handle
    % Class with constants data
    %
    % Attributes:
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     G (float): Standard gravity [m/s2]
    %     NA (float): Avogadro's number [molecule/mol]
    %     release (char): Release of the Combustion Toolbox
    %     date (char): Date of the release
    %
    % Examples:
    %     * R0 = Constants.R0
    %     * g = Constants.G
    %     * release = Constants.release
    
    properties (Constant)
        R0      = 8.31446261815324 % Universal gas constant [J/(K mol)]
        G       = 9.80665          % Standard gravity [m/s2]
        NA      = 6.0221415e23     % Avogadro's number [molecule/mol]
        release = 'v1.1.0beta'     % Release version
        date    = '03 Aug 2024'    % Release date
    end
    
end
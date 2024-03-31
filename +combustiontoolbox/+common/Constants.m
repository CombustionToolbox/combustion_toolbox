classdef Constants < handle
    % Initialize struct with constants data
    %
    % Attributes:
    %     release (char): Release of the Combustion Toolbox
    %     date (char): Date of the release
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     G (float): Standard gravity [m/s2]
    %
    % Returns:
    %     self (struct): Struct with constants data
    
    properties (Constant)
        R0      = 8.31446261815324 % Universal gas constant [J/(K mol)]
        G       = 9.80665          % Standard gravity [m/s2]
        NA      = 6.0221415e23     % Avogadro's number [molecule/mol]
    end
    
end
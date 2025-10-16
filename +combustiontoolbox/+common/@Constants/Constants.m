classdef Constants < handle
    % Class with constants data
    %
    % Attributes:
    %     R0 (float): Universal gas constant [J/(K mol)]
    %     G (float): Standard gravity [m/s2]
    %     NA (float): Avogadro's number [molecule/mol]
    %     KB (float): Boltzmann constant [J/K]
    %     HBAR (float): Planck constant [J-s]
    %     E (float): Electron charge [C]
    %     E0 (float): Vacuum permittivity [F/m]
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
        KB      = 1.38064852e-23   % Boltzmann constant [J/K]
        HBAR    = 6.626070041e-34  % Planck constant [J-s]
        E       = 1.602176634e-19  % Electron charge [C]
        E0      = 8.854187818e-12  % Vacuum permittivity [F/m]
        ME      = 9.10938356e-31   % Electron mass [kg]
        release = 'v1.2.8'         % Release version
        date    = '16 Oct 2025'    % Release date
    end
    
end

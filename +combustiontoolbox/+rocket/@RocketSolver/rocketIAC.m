function [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, varargin)
    % Routine that computes the propellant rocket performance using the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (RocketSolver): RocketSolver object
    %     mix1 (Mixture): Properties of the initial mixture
    %
    % Optional Args:
    %     * mix2_c (Mixture): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     * mix3 (Mixture): Properties of the mixture at the throat (previous calculation)
    %     * mix4 (Mixture): Properties of the mixture at the given exit points (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the initial mixture
    %     * mix2_c (Mixture): Properties of the mixture at the outlet of the chamber
    %     * mix3 (Mixture): Properties of the mixture at the throat
    %     * mix4 (Mixture): Array of mixtures objects at the given exit points
    %
    % Examples:
    %     * [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1)
    %     * [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1)
    %     * [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, mix2_c)
    %     * [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, mix2_c, mix3)
    %     * [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, mix2_c, mix3, mix4)
    

    % Unpack additional input parameters
    if nargin > 2
        mix2_c = copy(varargin{1});
        setProperties(mix2_c, 'temperature', mix1.T, 'pressure', mix1.p);
    else
        mix2_c = copy(mix1);
    end

    if nargin > 3
        mix3 = copy(varargin{2});
    else
        mix3 = [];
    end

    if nargin > 4
        mix4 = copy(varargin{3});
    else
        mix4 = [];
    end

    % Compute chemical equilibria at different stations of the rocket
    mix2_c = computeChamberIAC(obj, mix2_c);
    mix3 = computeThroatIAC(obj, mix2_c, mix3);
    mix4 = computeExit(obj, mix2_c, mix3, mix4, mix1.areaRatio);
    
    % Initial velocity of the gas
    mix1.u = 0; mix1.uShock = 0;

    % Velocity at the outlet of the chamber
    mix2_c.u = 0; mix2_c.uShock = 0;

    % Compute rocket parameters
    [mix3, mix4] = rocketParameters(obj, mix2_c, mix3, mix4);
end
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
    %     * FLAG_PERFORMANCE (bool): Flag to compute rocket parameters (default = true)
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
    
    % Initialization
    mix2_c = copy(mix1);
    FLAG_PERFORMANCE = true;

    % Unpack additional input parameters
    switch nargin
        case 3
            mix2_c_guess = varargin{1};
        case 4
            mix2_c_guess = varargin{1};
            mix3_guess = varargin{2};
        case 5
            mix2_c_guess = varargin{1};
            mix3_guess = varargin{2};
            mix4_guess = varargin{3};
        case 6
            mix2_c_guess = varargin{1};
            mix3_guess = varargin{2};
            mix4_guess = varargin{3};
            FLAG_PERFORMANCE = varargin{4};
        otherwise
            mix2_c_guess = [];
            mix3_guess = [];
            mix4_guess = [];
    end

    % Compute chemical equilibria at different stations of the rocket
    mix2_c = rocketChamberIAC(obj, mix2_c, mix2_c_guess);
    mix3 = rocketThroatIAC(obj, mix2_c, mix3_guess);
    mix4 = rocketExit(obj, mix2_c, mix3, mix4_guess, mix1.areaRatio);
    
    % Initial velocity of the gas
    mix1.u = 0; mix1.uShock = 0; mix1.mach = 0;
    mix1.I_sp = 0; mix1.I_vac = 0;
        
    % Velocity at the outlet of the chamber
    mix2_c.u = 0; mix2_c.uShock = 0; mix2_c.mach = 0;
    mix2_c.I_sp = 0; mix2_c.I_vac = 0;

    % Compute rocket parameters
    if ~FLAG_PERFORMANCE
        return
    end

    [mix3, mix4] = obj.rocketParameters(mix2_c, mix3, mix4);
end
function [mix1, mix2_inj, mix2_c, mix3, mix4] = rocket_performance(self, mix1, Aratio, varargin)
    % Routine that computes the propellant rocket performance
    %
    % Methods implemented:
    %   * Infinite-Area-Chamber (IAC) 
    %   * Finite-Area-Chamber (FAC) - NOT YET
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     Aratio (struct): Ratio area_exit / area_throat
    %
    % Optional Args:
    %     - mix2 (struct): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     - mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct): Properties of the initial mixture
    %     - mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     - mix3 (struct): Properties of the mixture at the throat
    %     - mix4 (struct): Properties of the mixture at the given exit points

    % Assign values
    if nargin > 3, mix2_inj = get_struct(varargin{1}); else, mix2_inj = []; end
    if nargin > 4, mix2_c = get_struct(varargin{2}); else, mix2_c = []; end
    if nargin > 5, mix3 = get_struct(varargin{3}); else, mix3 = []; end
    if nargin > 6, mix4 = get_struct(varargin{4}); else, mix4 = []; end
    % Compute chemical equilibria at different points of the rocket
    % depending of the model selected
    [mix2_inj, mix2_c, mix3, mix4] = solve_model_rocket(self, mix1, mix2_inj, mix2_c, mix3, mix4, Aratio);
    % Initial velocity of the gas
    mix1.u = 0; mix1.v_shock = 0;
    % Compute rocket parameters
    if self.PD.FLAG_IAC
        % Velocity at the outlet of the chamber
        mix2_c.u = 0; mix2_c.v_shock = 0;
        [mix3, mix4] = compute_rocket_parameters(mix2_c, mix3, self.C.gravity, mix4);
    else
        % Velocity at the injector
        mix2_inj.u = 0; mix2_inj.v_shock = 0;
        [mix3, mix2_c, mix4] = compute_rocket_parameters(mix2_inj, mix3, self.C.gravity, mix2_c, mix4);
    end
end

% SUB-PASS FUNCTIONS
function str = get_struct(var)
    try
        str = var{1,1};
    catch
        str = var;
    end
end
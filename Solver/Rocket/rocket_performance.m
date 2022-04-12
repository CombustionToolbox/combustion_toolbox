function [mix1, mix2, mix3] = rocket_performance(self, mix1, Aratio, varargin)
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
    if nargin > 3, mix2 = get_struct(varargin{1}); else, mix2 = []; end
    if nargin > 4, mix3 = get_struct(varargin{2}); else, mix3 = []; end
    if nargin > 5, mix4 = get_struct(varargin{3}); else, mix4 = []; end
    % Compute chemical equilibria at different points of the rocket
    % depending of the model selected
    [mix2_1, mix2, mix3, mix4] = solve_model(self, mix1, mix2, mix3, mix4, Aratio);
    % Velocity at the inlet and outlet of the chamber
    mix1.u = 0; mix1.v_shock = 0;
    mix2.u = 0; mix2.v_shock = 0;
    % Compute rocket parameters
    [mix3, mix4] = compute_rocket_parameters(mix2, mix3, self.C.gravity, mix4);
end

% SUB-PASS FUNCTIONS
function str = get_struct(var)
    try
        str = var{1,1};
    catch
        str = var;
    end
end

function [mix2_1, mix2, mix3, mix4] = solve_model(self, mix1, mix2, mix3, mix4, Aratio)
    % Compute chemical equilibria at different points of the rocket
    % depending of the model selected

    if self.PD.FLAG_IAC
        mix2_1 = [];
        mix2 = compute_chamber_IAC(self, mix1, mix2);
        mix3 = compute_throat_IAC(self, mix2, mix3);
        mix4 = compute_exit_IAC(self, mix2, mix3, mix4, Aratio);
    else
        [mix2_1, mix2] = compute_chamber_FAC(self, mix1, mix2);
        mix3 = compute_throat_FAC(self, mix2, mix3);
        mix4 = compute_exit_FAC(self, mix2, mix3, mix4, Aratio);
    end
end
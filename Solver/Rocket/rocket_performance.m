function [mix1, mix2, mix3] = rocket_performance(self, mix1, varargin)
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
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     mix1 (struct): Properties of the initial mixture
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     mix3 (struct): Properties of the mixture at the throat

    % Assign values
    if nargin == 3, mix2 = get_struct(varargin{1}); else, mix2 = []; end
    if nargin == 4, mix3 = get_struct(varargin{2}); else, mix3 = []; end
    % Compute chemical equilibria at the exit of the chamber (HP)
    mix2 = compute_chamber(self, mix1, mix2);
    % Compute chemical equilibria at throat (SP)
    mix3 = compute_throat(self, mix2, mix3);
    % Velocity at the inlet and outlet of the chamber
    mix1.u = 0; mix1.v_shock = 0; 
    mix2.u = 0; mix2.v_shock = 0;
    % Compute rocket parameters
    mix3 = compute_rocket_parameters(mix2, mix3);
end

% SUB-PASS FUNCTIONS
function str = get_struct(var)
    try
        str = var{1,1};
    catch
        str = var;
    end
end

function mix2 = compute_chamber(self, mix1, mix2)
    % Compute chemical equilibria at the exit of the chamber (HP)
    self.PD.ProblemType = 'HP';
    mix2 = compute_chemical_equilibria(self, mix1, mix1.p, mix2);
end

function mix3 = compute_throat(self, mix2, mix3)
    % Compute chemical equilibria at the throat (SP)
    self = set_prop(self, 'TR', mix2.T, 'pR', mix2.p);
    self.PD.S_Fuel = self.S.LS;
    self.PD.N_Fuel = moles(mix2)';
    self.PD.ProblemType = 'SP';
    mix3 = compute_IAC_model(self, mix2, mix3);
end
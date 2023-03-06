function [mix2_inj, mix2_c, mix3, mix4] = solve_model_rocket(self, mix1, mix2_inj, mix2_c, mix3, mix4, Aratio)
    % Compute chemical equilibria at different points of the rocket
    % depending of the model selected
    %
    % Methods implemented:
    %   * Infinite-Area-Chamber (IAC)
    %   * Finite-Area-Chamber (FAC)
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     mix2_inj (struct): Properties of the mixture at the injector [only FAC] (previous calculation)
    %     mix2_c (struct): Properties of the mixture at the outlet of the chamber (previous calculation)
    %     mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %     mix4 (struct): Properties of the mixture at the given exit points (previous calculation)
    %     Aratio (float): Area ratios [-]
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix2_1 (struct): Properties of the mixture at injector of the chamber (only FAC)
    %     * mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     * mix3 (struct): Properties of the mixture at the throat
    %     * mix4 (struct): Properties of the mixture at the given exit points

    if self.PD.FLAG_IAC
        mix2_inj = [];
        mix2_c = compute_chamber_IAC(self, mix1, mix2_c);
        mix3 = compute_throat_IAC(self, mix2_c, mix3);
        mix4 = compute_exit(self, mix2_c, mix3, mix4, Aratio);
    else
        [mix2_inj, mix2_c, mix3] = compute_FAC(self, mix1, mix2_inj, mix2_c, mix3);
        mix4 = compute_exit(self, mix2_c, mix3, mix4, Aratio, mix2_inj);
    end

end

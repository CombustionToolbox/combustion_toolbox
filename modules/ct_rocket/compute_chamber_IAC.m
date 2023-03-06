function mix2 = compute_chamber_IAC(self, mix1, mix2)
    % Compute chemical equilibria at the exit of the chamber (HP) using
    % the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber (previous calculation)
    %
    % Returns:
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber

    % Definitions
    self.PD.ProblemType = 'HP';
    % Compute chemical equilibria at the exit of the chamber (HP)
    mix2 = equilibrate(self, mix1, mix1.p, mix2);
    % Set A_chamber/A_throat
    mix2.Aratio = Inf;
end

function mix2 = compute_chamber_IAC(self, mix1, mix2)
    % Compute chemical equilibria at the exit of the chamber (HP) using
    % the Infinite-Area-Chamber (IAC) model
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix2 (struct): Properties of the mixture at the outlet of the chamber
    %     mix3 (struct): Properties of the mixture at the throat (previous calculation)
    %
    % Returns:
    %     mix3 (struct): Properties of the mixture at the throat
    
    % Definitions
    self.PD.ProblemType = 'HP';
    % Compute chemical equilibria at the exit of the chamber (HP)
    mix2 = compute_chemical_equilibria(self, mix1, mix1.p, mix2);
end
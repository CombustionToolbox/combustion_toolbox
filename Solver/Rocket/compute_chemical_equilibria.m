function mix2 = compute_chemical_equilibria(self, mix1, pP, mix2)
    % Compute chemical equilibria for the 2 given thermodynamic states,
    % e.g., enthalpy-pressure (HP)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float):    Pressure [bar]
    %     mix2 (struct): Properties of the final mixture (previous calculation)
    % 
    % Returns:
    %     mix2 (struct): Properties of the final mixture

    if isempty(mix2)
        mix2 = equilibrate(self, mix1, pP);
    else
        mix2 = equilibrate(self, mix1, pP, mix2);
    end
end

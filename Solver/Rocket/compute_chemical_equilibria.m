function mix2 = compute_chemical_equilibria(self, mix1, pP, mix2)
    % Compute chemical equilibria for the 2 given thermodynamic states,
    % e.g., enthalpy-pressure (HP)
    if isempty(mix2)
        mix2 = equilibrate(self, mix1, pP);
    else
        mix2 = equilibrate(self, mix1, pP, mix2);
    end
end

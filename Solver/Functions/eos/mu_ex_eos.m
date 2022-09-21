function chemical_potential_ex = mu_ex_eos(self, T, p, v_molar, a, b)
    % Compute non ideal contribution (excess) of the chemical potential
    % assuming ideal Equation of State [J/mol]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): number of moles of the mixture in gaseous phase [mol]
    %     temperature (float): temperature of the mixture [K]
    %     volume (float): volume of the mixture [m3]
    % 
    % Returns:
    %     pressure (float): pressure of the mixture [Pa]
    
    % Constants
    R0 = self.C.R0;
    % Change of units
    p = convert_bar_to_Pa(p);
    % Definitions
    Z = p * v_molar / (R0 * T);
    % Compute excess of the chemical potential
    chemical_potential_ex = p * v_molar - R0 * T * (1 + log(Z - p * b / (R0 * T))) - a / (2 * sqrt(2) * b) * log((v_molar + (1 + sqrt(2)) * b) / (v_molar + (1 - sqrt(2)) * b));
end
function chemical_potential_ex = mu_ex_ideal(self, moles, temperature, volume)
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

    chemical_potential_ex = 0;
end
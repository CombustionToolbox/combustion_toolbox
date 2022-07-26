function chemical_potential_imp = mu_imp_ideal(self, moles, temperature, volume)
    % Compute non ideal contribution of the chemical potential assuming 
    % ideal Equation of States [J/mol]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): number of moles of the mixture in gaseous phase [mol]
    %     temperature (float): temperature of the mixture [K]
    %     volume (float): volume of the mixture [m3]
    % 
    % Returns:
    %     pressure (float): pressure of the mixture [Pa]

    chemical_potential_imp = 0;
end
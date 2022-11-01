function chemical_potential_ex = mu_ex_vanderwaals(self, moles, temperature, volume)
    % Compute non ideal contribution (excess) of the chemical potential
    % assuming Van der Waal's Equation of State [J/mol]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): number of moles of the mixture in gaseous phase [mol]
    %     temperature (float): temperature of the mixture [K]
    %     volume (float): volume of the mixture [m3]
    % 
    % Returns:
    %     pressure (float): pressure of the mixture [Pa]
    
    % Definitions
    Nmoles = length(moles);
    molesGas = sum(moles);
    % Compute coefficients
    a = compute_cofficient(moles, molesGas, Nmoles, self.PD.EOS.a);
    b = compute_cofficient(moles, molesGas, Nmoles, self.PD.EOS.b);
    
    volume_molar = volume * self.C.NA;
    pressure = self.C.R0 * temperature / (volume_molar - b) - a / volume_molar^2;
    chemical_potential_ex = 0;
end

function value = compute_cofficient(moles, molesGas, coefficients, Nmoles)
    % Comptue coefficient a
    value = 0;
    for i = 1:Nmoles
        for j = 1:Nmoles
            value = value + moles(i) * moles(j) * sqrt(coefficients(i) * coefficients(j));
        end
    end
    value = value / molesGas^2;
end
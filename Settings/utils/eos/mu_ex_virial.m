function chemical_potential_ex = mu_ex_virial(self, moles, temperature, volume)
    % Compute non ideal contribution (excess) of the chemical potential
    % assuming Virial Equation of State [J/mol]
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
    a_ij = self.PD.EOS.a;
    b_j = self.PD.EOS.b;
    % Compute coefficients
    a = compute_cofficient_pair(moles, molesGas, self.PD.EOS.a, Nmoles);
%     b = compute_cofficient_single(moles, molesGas, );
    c = sum(moles .* a_ij - b_j);
    chemical_potential_ex = self.C.R0 * temperature * volume^2 * (-a^2 - 2*b + 4*a / molesGas * c);
end

function value = compute_cofficient_pair(moles, molesGas, coefficients, Nmoles)
    % Comptue coefficient a
    value = 0;
    for i = 1:Nmoles
        for j = 1:Nmoles
            value = value + moles(i) * moles(j) * sqrt(coefficients(i) * coefficients(j));
        end
    end
    value = value / molesGas^2;
end

function value = compute_cofficient_single(moles, molesGas, coefficients)
    % Comptue coefficient a
    value = dot(moles, coefficients) / molesGas;
end
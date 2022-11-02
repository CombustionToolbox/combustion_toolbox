function chemical_potential_ex = mu_ex_eos(self, Xi, T, p, V, a_mix, b_mix, a, b)
    % Compute non ideal contribution (excess) of the chemical potential
    % assuming cubic Equation of State (EoS) [J/mol]
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     T (float): Temperature of the mixture [K]
    %     p (float): Pressure of the mixture [bar]
    %     V (float): Molar volume of the mixture [m3/mol]
    %     a_mix (float): Atraction factor mixture of the cubic EoS 
    %     b_mix (float): Repulsion factor mixture of the cubic EoS
    %     a (float): Atraction factor components of the cubic EoS 
    %     b (float): Repulsion factor components of the cubic EoS
    %
    % Returns:
    %     chemical_potential_ex (float): chemical potential excess [J/mol]
    
    % Constants
    R0 = self.C.R0;
    % Change of units
    p = convert_bar_to_Pa(p);
    % Definitions
    Z = p * V / (R0 * T);
    % Compute excess of the chemical potential
    chemical_potential_ex(:, 1) = p * V - R0 * T * (1 + log(Z - p * b_mix / (R0 * T))) - a_mix / (2 * sqrt(2) * b_mix) * (R0 * T * (2 * a * Xi/ a_mix) - b / b_mix) * log((V + (1 + sqrt(2)) * b_mix) / (V + (1 - sqrt(2)) * b_mix));
end
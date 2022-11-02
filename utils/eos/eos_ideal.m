function V = eos_ideal(self, T, p, varargin)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     T (float): temperature of the mixture [K]
    %     p (float): pressure of the mixture [Pa]
    % 
    % Returns:
    %     V (float): molar volume of the mixture [m3/mol]

    V = self.C.R0 * T ./ p;
end
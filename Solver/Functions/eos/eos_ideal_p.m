function p = eos_ideal_p(self, n, T, v, varargin)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     n (float): number of moles of the mixture in gaseous phase [mol]
    %     T (float): temperature of the mixture [K]
    %     v (float): volume of the mixture [m3]
    % 
    % Returns:
    %     p (float): pressure of the mixture [Pa]

    p = (n .* self.C.R0 .* T) ./ v;
end
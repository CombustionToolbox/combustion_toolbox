function p = eos_ideal_p(self, n, T, v, varargin)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     n (float): Number of moles of the mixture in gaseous phase [mol]
    %     T (float): Temperature of the mixture [K]
    %     v (float): Volume of the mixture [m3]
    % 
    % Returns:
    %     p (float): Pressure of the mixture [Pa]

    p = (n .* self.C.R0 .* T) ./ v;
end
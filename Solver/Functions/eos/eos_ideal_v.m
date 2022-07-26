function volume = eos_ideal_v(self, moles, temperature, pressure)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): number of moles of the mixture in gaseous phase [mol]
    %     temperature (float): temperature of the mixture [K]
    %     pressure (float): pressure of the mixture [Pa]
    % 
    % Returns:
    %     volume (float): volume of the mixture [m3]

    volume = (moles .* self.C.R0 .* temperature) ./ pressure;
end
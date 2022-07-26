function pressure = eos_ideal(self, moles, temperature, volume)
    % Compute pressure considering ideal Equation of State (EoS)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): number of moles of the mixture in gaseous phase [mol]
    %     temperature (float): temperature of the mixture [K]
    %     volume (float): volume of the mixture [m3]
    % 
    % Returns:
    %     pressure (float): pressure of the mixture [Pa]

    pressure = (moles .* self.C.R0 .* temperature) ./ volume;
end
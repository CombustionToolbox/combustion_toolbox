function value = get_transformation(self, field)
    % Get the corresponding value of the field in Problem Description (PD)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     field (str):   Fieldname in Problem Description (PD)
    %
    % Returns:
    %     value (float): Value/s assigned to the field
    
    value = self.PD.(field).value;
end
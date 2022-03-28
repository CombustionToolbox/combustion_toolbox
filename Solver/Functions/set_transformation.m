function self = set_transformation(self, field, value)
    % Set the corresponding value of the field in Problem Description (PD)
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     field (str):   Fieldname in Problem Description (PD)
    %     value (float): Value/s to assign to the field
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    self.PD.(field).value = value;
end
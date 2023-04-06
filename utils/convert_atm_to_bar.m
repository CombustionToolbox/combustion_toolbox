function value = convert_atm_to_bar(value)
    % Convert pressure in [atm] units to [bar]
    %
    % Args:
    %     value (float): Pressure value(s) in [atm]
    %
    % Returns:
    %     value (float): Pressure value(s) in [bar]
    %
    % Example:
    %     value = convert_atm_to_bar(1)

    value = value * 1.01325;
end

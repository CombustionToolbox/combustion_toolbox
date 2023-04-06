function value = convert_bar_to_atm(value)
    % Convert pressure in [bar] units to [atm]
    %
    % Args:
    %     value (float): Pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): Pressure value(s) in [atm]
    %
    % Example:
    %     value = convert_bar_to_atm(1.01325)

    value = value / 1.01325;
end

function value = convert_bar_to_atm(value)
    % Convert pressure in [bar] units to [atm]
    %
    % Args:
    %     value (float): pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): pressure value(s) in [atm]

    value = value / 1.01325;
end

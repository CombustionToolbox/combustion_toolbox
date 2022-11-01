function value = convert_atm_to_bar(value)
    % Convert pressure in [atm] units to [bar]
    %
    % Args:
    %     value (float): pressure value(s) in [atm]
    %
    % Returns:
    %     value (float): pressure value(s) in [bar]

    value = value * 1.01325;
end

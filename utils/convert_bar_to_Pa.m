function value = convert_bar_to_Pa(value)
    % Convert pressure in [bar] units to [Pa]
    %
    % Args:
    %     value (float): pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): pressure value(s) in [Pa]

    value = value * 1e5;
end

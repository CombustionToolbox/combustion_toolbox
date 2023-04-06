function value = convert_bar_to_Pa(value)
    % Convert pressure in [bar] units to [Pa]
    %
    % Args:
    %     value (float): Pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): Pressure value(s) in [Pa]
    %
    % Example:
    %     value = convert_bar_to_Pa(1)

    value = value * 1e5;
end

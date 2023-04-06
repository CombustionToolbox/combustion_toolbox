function value = convert_Pa_to_bar(value)
    % Convert pressure in [Pa] units to [bar]
    %
    % Args:
    %     value (float): Pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): Pressure value(s) in [bar]
    %
    % Example:
    %     value = convert_Pa_to_bar(1e5)

    value = value * 1e-5;
end

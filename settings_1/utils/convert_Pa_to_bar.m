function value = convert_Pa_to_bar(value)
    % Convert pressure in [Pa] units to [bar]
    %
    % Args:
    %     value (float): pressure value(s) in [bar]
    %
    % Returns:
    %     value (float): pressure value(s) in [bar]

    value = value * 1e-5;
end

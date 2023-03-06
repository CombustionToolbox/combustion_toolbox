function point = get_point(x_vector, f_vector)
    % Get point using the regula falsi method
    %
    % Args:
    %     x_vector (float):  Guess temperature [K]
    %     f_vector (struct): evaluated functions [kJ] (HP, EV) or [kJ/K] (SP, SV)
    % Returns:
    %     point (float): Point of the function [K]

    point = (f_vector(2) * x_vector(1) - f_vector(1) * x_vector(2)) / (f_vector(2) - f_vector(1));
end
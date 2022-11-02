function point = get_point(x_vector, g_vector)
    % Get point of the fixed point function
    %
    % Args:
    %     x_vector (float):  Guess temperature [K]
    %     g_vector (struct): Fixed points of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)
    % Returns:
    %     point (float): Point of the function [K]

    point =  x_vector(2) - (x_vector(2) - x_vector(1)) / (g_vector(2) - g_vector(1)) * g_vector(2);
end
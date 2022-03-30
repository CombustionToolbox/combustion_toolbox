function point = get_point_aitken(x0, g_vector)
    % Get fixed point of a function based on the chemical transformation using the Aitken acceleration method
    %
    % Args:
    %     x0 (float):        Guess temperature [K]
    %     g_vector (struct): Fixed points of the function [kJ] (HP, EV) or [kJ/K] (SP, SV)
    %
    % Returns:
    %     point (float): Point of the function [K]
    
    point = x0 - (g_vector(1) - x0)^2 / (g_vector(2) - 2*g_vector(1) + x0);   
end
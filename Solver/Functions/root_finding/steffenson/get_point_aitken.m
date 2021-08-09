function point = get_point_aitken(x, g_vector)
    point = x - (g_vector(1) - x)^2 / (g_vector(2) - 2*g_vector(1) + x);   
end
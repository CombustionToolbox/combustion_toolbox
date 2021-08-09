function point = get_point(x_vector, g_vector)
    point =  x_vector(2) - (x_vector(2) - x_vector(1)) / (g_vector(2) - g_vector(1)) * g_vector(2);
end
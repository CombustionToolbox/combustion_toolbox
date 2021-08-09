function x = steff_guess(self, strR, pP, attr_name)
    x_l = strR.T + 50;
    x_r = 1500;
            
    g_l = get_gpoint(self, strR, pP, attr_name, x_l);
    g_r = get_gpoint(self, strR, pP, attr_name, x_r);
    
    if g_l * g_r > 0 && g_l < g_r
        x = x_l - 50;
    elseif g_l * g_r > 0 || (np.isnan(g_l) && isnan(g_r))
        x = 2000;
    elseif isnan(g_l) && ~isnan(g_r)
        x = x_r - 100;
    elseif ~isnan(g_l) && isnan(g_r)
        x = x_l + 100;
    else
        x = get_point([x_l, x_r], [g_l, g_r]);
    end    
end
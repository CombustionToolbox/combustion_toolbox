function [x, ERR] = steff(self, strR, pP, attr_name, x)
    % Steffenson method for finding roots
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        ERR = 0;
        return
    end
    it = 0; g = 1.0; ERR = 1.0;
    
    while (abs(ERR) > self.TN.tol0 || abs(g) > self.TN.tol0) && it < self.TN.itMax
        it = it + 1;
        g = get_gpoint(self, strR, pP, attr_name, x);
        fx = abs(g - x);
        g_aux  = get_gpoint(self, strR, pP, attr_name, fx);
        fx2 = abs(g_aux - fx);
        if abs(fx2 - 2*fx + x) > self.TN.tol0
            x = get_point_aitken(x, [fx, fx2]);
        else
            x = fx;
        end

        ERR = max(abs((x - fx) / x), abs(g_aux - g));
    end
    print_error_root(it, self.TN.itMax, x, ERR);
end
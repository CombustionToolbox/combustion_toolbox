function [x, ERR] = newton(self, mix1, pP, attr_name, x0)
    % 1D Newton-Raphson method for finding roots
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        ERR = 0;
        return
    end
    it = 0; ERR = 1.0;
    while abs(ERR) > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        [f0, fprime0] = get_ratio_newton(self, mix1, pP, attr_name, x0);
        x = abs(x0 - f0 / fprime0);
        
        f = get_gpoint(self, mix1, pP, attr_name, x);
        ERR = max(abs((x - x0) / x), abs((f - f0) / f));
        x0 = x;
    end
    print_error_root(it, self.TN.itMax, x, ERR);
end

%%% SUB-PASS FUNCTIONS
function [f, fprime] = get_ratio_newton(self, mix1, pP, attr_name, x)
    mix2 = equilibrate_T(self, mix1, pP, x);
    f = mix2.(attr_name) - mix1.(attr_name);
    fprime = get_partial_derivative(self, mix2);
end
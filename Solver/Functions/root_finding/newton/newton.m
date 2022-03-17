function [x, ERR] = newton(self, mix1, pP, attr_name, x0)
    % 1D Newton-Raphson method for finding roots
    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        ERR = 0;
        return
    end
    it = 0; ERR = 1.0;
    while ERR > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        [f0, fprime0, frel] = get_ratio_newton(self, mix1, pP, attr_name, x0);
        x = abs(x0 - f0 / fprime0);
        
        ERR = max(abs((x - x0) / x), frel);
        x0 = x;
    end
    if ERR > self.TN.tol0
        fprintf('\n***********************************************************\n')
        fprintf('Newton method not converged\nCalling Steffensen-Aitken root finding algorithm\n')
        x0 = steff_guess(self, mix1, pP, attr_name);
        [x, ERR] = steff(self, mix1, pP, attr_name, x0);
    else
        print_error_root(it, self.TN.itMax, x, ERR);
    end
end

%%% SUB-PASS FUNCTIONS
function [f, fprime, frel] = get_ratio_newton(self, mix1, pP, attr_name, x)
    mix2 = equilibrate_T(self, mix1, pP, x);
    f = mix2.(attr_name) - mix1.(attr_name);
    fprime = get_partial_derivative(self, mix2);
    frel = abs(f / mix2.(attr_name));
end
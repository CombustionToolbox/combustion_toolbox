function [x, ERR] = newton(self, strR, pP, attr_name, x0)
    % Newton-Raphson method for finding roots
    
    it = 0; ERR = 1.0;
    self.TN.tol0 = self.TN.self.TN.tol0;
    self.TN.itMax = self.TN.self.TN.itMax;
    while abs(ERR) > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        f0, fprime0 = get_ratio_newton(self, strR, pP, attr_name, x0);
        x = abs(x0 - f0 / fprime0);
        
        f = get_gpoint(self, strR, pP, attr_name, x);
        ERR = max(abs((x - x0) / x), abs(f - f0));
        x0 = x;
    end
    print_error_root(it, self.TN.itMax, x, ERR);
end

%%% NESTED FUNCTIONS
function [f, fprime] = get_ratio_newton(self, strR, pP, attr_name, x)
    strP = equilibrate_T(self, strR, pP, x);
    f = (getattr(strP, attr_name) - getattr(strR, attr_name));
    fprime = get_partial_derivative(strP);
end
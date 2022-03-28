function [x, STOP] = newton(self, mix1, pP, field, x0)
    % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Newton-Raphson method
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float):    Pressure [bar]
    %     field (str):   Fieldname in Problem Description (PD)
    %     x0 (float):    Guess temperature [K]
    % Returns:
    %     x (float):     Temperature at equilibrium [K]
    %     STOP (float):  Relative error [-] 

    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        STOP = 0;
        return
    end

    it = 0; STOP = 1.0;

    while STOP > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        [f0, fprime0, frel] = get_ratio_newton(self, mix1, pP, field, x0);
        x = abs(x0 - f0 / fprime0);
        
        STOP = max(abs((x - x0) / x), frel);
        x0 = x;
    end
    if STOP > self.TN.tol0
        fprintf('\n***********************************************************\n')
        fprintf('Newton method not converged\nCalling Steffensen-Aitken root finding algorithm\n')
        x0 = steff_guess(self, mix1, pP, field);
        [x, STOP] = steff(self, mix1, pP, field, x0);
    else
        print_error_root(it, self.TN.itMax, x, STOP);
    end
end

%%% SUB-PASS FUNCTIONS
function [f, fprime, frel] = get_ratio_newton(self, mix1, pP, field, x)
    mix2 = equilibrate_T(self, mix1, pP, x);
    f = mix2.(field) - mix1.(field);
    fprime = get_partial_derivative(self, mix2);
    frel = abs(f / mix2.(field));
end
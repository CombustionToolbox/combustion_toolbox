function [x, STOP, guess_moles] = newton(self, mix1, pP, field, x0, guess_moles)
    % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the Newton-Raphson method
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float): Pressure [bar]
    %     field (str): Fieldname in Problem Description (PD)
    %     x0 (float): Guess temperature [K]
    %     guess_moles (float): Guess moles final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     - x (float): Temperature at equilibrium [K]
    %     - STOP (float): Relative error [-] 
    %     - guess_moles (struct): Guess moles final mixture

    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        STOP = 0;
        return
    end

    it = 0; STOP = 1.0;

    while STOP > self.TN.tol0 && it < self.TN.itMax
        it = it + 1;
        [f0, fprime0, frel, guess_moles] = get_ratio_newton(self, mix1, pP, field, x0, guess_moles);
        x = abs(x0 - f0 / fprime0);
        % Re-estimation of first derivative
%         [~, fprime0_2] = get_ratio_newton(self, mix1, pP, field, x, guess_moles);
        % Pseudo re-estimation of first derivative
        fprime0_2 = fprime0 * x0/x;
        x = abs(x0 - 2*f0 / (fprime0 + fprime0_2));

        STOP = max(abs((x - x0) / x), frel);
        x0 = x;
        % Debug       
%         aux_x(it) = x;
%         aux_STOP(it) = STOP;
    end
%     debug_plot_error(it, aux_STOP, aux_x);
    if STOP > self.TN.tol0
        fprintf('\n***********************************************************\n')
        fprintf('Newton method not converged\nCalling Steffensen-Aitken root finding algorithm\n')
        x0 = steff_guess(self, mix1, pP, field);
        [x, STOP] = steff(self, mix1, pP, field, x0, []);
    else
        print_error_root(it, self.TN.itMax, x, STOP);
    end
end

%%% SUB-PASS FUNCTIONS
function [f, fprime, frel, guess_moles] = get_ratio_newton(self, mix1, pP, field, x, guess_moles)
    mix2 = equilibrate_T(self, mix1, pP, x, guess_moles);
    f = mix2.(field) - mix1.(field);
    fprime = get_partial_derivative(self, mix2);
    frel = abs(f / mix2.(field));
    % Update guess moles
    guess_moles = mix2.N * mix2.Xi;
end
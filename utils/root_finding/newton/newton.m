function [x, STOP, guess_moles] = newton(self, mix1, pP, field, x0, guess_moles)
    % Find the temperature [K] (root) for the set chemical transformation at equilibrium using the second-order Newton-Raphson method
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
    %     * x (float): Temperature at equilibrium [K]
    %     * STOP (float): Relative error [-] 
    %     * guess_moles (struct): Guess moles final mixture

    if any(strcmpi(self.PD.ProblemType, {'TP', 'TV'}))
        x = get_transformation(self, 'TP');
        STOP = 0;
        return
    end
    
    % Initialization
    it = 0; STOP = 1.0;

    % Loop
    while STOP > self.TN.tol0 && it < self.TN.itMax
        % Update iteration number
        it = it + 1;
        % Get the residual of f, its derivative with temperature, and the
        % relative value of the residual
        [f0, fprime0, frel, guess_moles] = get_ratio_newton(self, mix1, pP, field, x0, guess_moles);
        % Compute solution
        x = abs(x0 - f0 / fprime0);
        % Compute stop criteria
        STOP = max(abs((x - x0) / x), frel);
        % Update solution
        x0 = x;
        % Debug       
        % aux_x(it) = x;
        % aux_STOP(it) = STOP;
    end

    % debug_plot_error(it, aux_STOP);
    if STOP > self.TN.tol0
        fprintf('\n***********************************************************\n')
        fprintf('Newton method not converged\nCalling Newton-Steffensen root finding algorithm\n')
        x0 = regula_guess(self, mix1, pP, field);
        [x, STOP] = nsteff(self, mix1, pP, field, x0, []);
    else
        print_error_root(it, self.TN.itMax, x, STOP);
    end
    
end

% SUB-PASS FUNCTIONS
function [f, fprime, frel, guess_moles] = get_ratio_newton(self, mix1, pP, field, x, guess_moles)
    % Get the residual of f, its derivative with temperature, and the
    % relative value of the residual
    try
        mix2 = equilibrate_T(self, mix1, pP, x, guess_moles);
    catch
        mix2 = equilibrate_T(self, mix1, pP, x);
    end

    % Calculate residual of f = 0
    f = mix2.(field) - mix1.(field);
    % Calculate partial derivative of f with temperature
    fprime = get_partial_derivative(self, mix2);
    % Get relative value of the residual
    frel = abs(f / mix2.(field));
    % Update guess moles
    guess_moles = mix2.N * mix2.Xi;
end
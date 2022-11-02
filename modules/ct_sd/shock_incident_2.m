function [mix1, mix2] = shock_incident_2(self, mix1, u1, varargin)
    % Compute pre-shock and post-shock states of a planar incident shock wave
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     u1 (float):    Pre-shock velocity [m/s]
    %
    % Optional Args:
    %     mix2 (struct): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     - mix1 (struct): Properties of the mixture in the pre-shock state
    %     - mix2 (struct): Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2, guess_moles] = unpack(self, mix1, u1, varargin);
    % Abbreviations
    TN = self.TN;
    % Miscellaneous
    self.PD.ProblemType = 'TP';
    % Definitions
    FLAG_FAST = TN.FLAG_FAST;
    lambda = 0.5; % Correction factor
    T1 = temperature(mix1); % [K]
    h1 = enthalpy_mass(mix1) * 1e3; % [J/kg]
    Rg1 = self.C.R0 / meanMolecularWeight(mix1) * 1e3; % [J/K];
    % Solve shock incident
    [mix2, STOP] = solve_shock_incident(mix2, FLAG_FAST);
    % If error, repeat without guess
    if isnan(STOP)
        fprintf('Recalculating: %.2f [m/s]\n', u1);
        guess_moles = [];
        FLAG_FAST = false;
        [mix2, STOP] = solve_shock_incident(mix2, FLAG_FAST);
    end

    % Check convergence
    print_convergence(STOP, TN.tol_shocks, temperature(mix2));
    % Save state
    mix2 = save_state(mix1, mix2, STOP);

    % NESTED-FUNCTIONS
    function [mix2, STOP] = solve_shock_incident(mix2, FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = TN.it_shocks; STOP = 2;
        % Initial estimate of post-shock state
        [mix2, u20] = get_guess(self, mix1, mix2);
        % Get properties
        [cp2, h2, Rg2, rho2] = get_properties(self, mix2);
        % Check FLAG
        if ~FLAG_FAST, guess_moles = []; end
        % Loop
        u2 = u20;

        while STOP > TN.tol_shocks && it < itMax
            it = it + 1;
            % Obtain f0 and fprime0
            f0 = h2 - h1 + 0.5 * (u2^2 - u1^2);
            fprime0 = cp2 / Rg2 * (Rg1 * T1 / u1 + u1 - 2 * u2) + u2;
            frel = abs(f0 / u2);
            % Apply correction
            u2 = u20 - lambda * f0 / fprime0; % [m/s]
            % Update post-shock state (temperature, pressure, mixture)
            T2 = u2 / Rg2 * (Rg1 * T1 / u1 + u1 - u2); % [K]
            p2 = convert_Pa_to_bar(rho2 * Rg2 * T2); % [bar]
            mix2 = equilibrate_T(self, mix1, p2, T2, guess_moles);
            % Update properties
            [cp2, h2, Rg2, rho2] = get_properties(self, mix2);
            % Compute STOP criteria
            STOP = compute_STOP(u2, u20, frel);
            % Update post-shock velocity u20
            u20 = u2;
            % Debug
            %             aux_lambda(it) = lambda;
            %             aux_STOP(it) = STOP;
        end

        %         debug_plot_error(it, aux_STOP, aux_lambda);
    end

end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2, guess_moles] = unpack(self, mix1, u1, x)
    % Unpack input data
    mix1.u = u1; % velocity pre-shock [m/s] - laboratory fixed
    mix1.v_shock = u1; % velocity pre-shock [m/s] - shock fixed

    try
        mix2 = x{1};
        guess_moles = mix2.Xi * mix2.N;
    catch
        mix2 = [];
        guess_moles = [];
    end

end

function [mix2, u2] = get_guess(self, mix1, mix2)

    if isempty(mix2)
        % Estimate post-shock state considering h2 = h1 + u1^2 / 2
        M1 = mix1.u / mix1.sound;
        p2p1 = (2 * mix1.gamma * M1^2 - mix1.gamma + 1) / (mix1.gamma + 1);
        mix1.h = (enthalpy_mass(mix1) * 1e3 + velocity_relative(mix1)^2/2) * mass(mix1) * 1e-3; % [kJ]
        self.PD.ProblemType = 'HP';
        mix2 = equilibrate(self, mix1, p2p1 * mix1.p);
    end

    u2 = mix1.rho * mix1.u / mix2.rho;
end

function [cp, h, Rg, rho] = get_properties(self, mix)
    cp = cp_mass(mix) * 1e3; % [J/kg-K]
    h = enthalpy_mass(mix) * 1e3; % [J/kg]
    Rg = self.C.R0 / meanMolecularWeight(mix) * 1e3; % [J/K];
    rho = density(mix); % [kg/m3]
end

function mix2 = save_state(mix1, mix2, STOP)
    mix2.v_shock = mix1.u * mix1.rho / mix2.rho;
    mix2.u = mix1.u - mix2.v_shock; % velocity postshock [m/s] - laboratory fixed
    mix2.error_problem = STOP;
end

function STOP = compute_STOP(x, x0, frel)
    % compute stop condition
    STOP = max(abs((x - x0) / x), frel);
end

function print_convergence(STOP, TOL, T)

    if STOP > TOL
        fprintf('***********************************************************\n')
        fprintf('Convergence error: %.2f\n', STOP);
    end

    if T > 2e4
        fprintf('***********************************************************\n')
        fprintf('Validity of the next results compromise\n')
        fprintf('Thermodynamic properties fitted to 20000 K\n');
    end

end

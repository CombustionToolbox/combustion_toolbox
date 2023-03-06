function [mix1, mix2] = shock_incident(self, mix1, u1, varargin)
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
    %     * mix1 (struct): Properties of the mixture in the pre-shock state
    %     * mix2 (struct): Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2, guess_moles] = unpack(self, mix1, u1, varargin);

    if self.TN.FLAG_TCHEM_FROZEN
        STOP = 0;
        [~, p2p1, T2T1, ~, ~] = shock_ideal_gas(mix1.gamma, u1 / mix1.sound);
        T2 = T2T1 * mix1.T;
        p2 = p2p1 * convert_bar_to_Pa(mix1.p);
        mix2 = save_state(self, mix1, T2, p2, STOP);
        return
    end
    % Abbreviations
    C = self.C;
    TN = self.TN;
    % Constants
    R0 = C.R0; % Universal gas constant [J/(mol-K)]
    % Definitions
    FLAG_FAST = TN.FLAG_FAST;
    % Solve shock incident
    [T2, p2, STOP] = solve_shock_incident(FLAG_FAST);
    if STOP > TN.tol_shocks
        % If solution has not converged, repeat without composition estimate
        fprintf('Recalculating: %.2f [m/s]\n', u1);
        guess_moles = [];
        FLAG_FAST = false;
        [T2, p2, STOP] = solve_shock_incident(FLAG_FAST);
    end

    % Check convergence
    print_convergence(STOP, TN.tol_shocks, T2);
    % Save state
    mix2 = save_state(self, mix1, T2, p2, STOP);

    % NESTED-FUNCTIONS
    function [T2, p2, STOP] = solve_shock_incident(FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = TN.it_shocks; STOP = 1.0;
        % Initial estimates of p2/p1 and T2/T1
        [p2, T2, p2p1, T2T1] = get_guess(self, mix1, mix2, TN);
        % Check FLAG
        if ~FLAG_FAST, guess_moles = []; end
        % Loop
        while STOP > TN.tol_shocks && it < itMax
            it = it + 1;
            % Construction of the Jacobian matrix and vector b
            [J, b, guess_moles] = update_system(self, mix1, p2, T2, R0, guess_moles, FLAG_FAST);
            % Solve of the linear system J*x = b
            x = J \ b;
            % Calculate correction factor
            lambda = relax_factor(x);
            % Apply correction
            [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda);
            % Apply antilog
            [p2, T2] = apply_antilog(mix1, log_p2p1, log_T2T1); % [Pa] and [K]
            % Update ratios
            p2p1 = p2 / (convert_bar_to_Pa(mix1.p));
            T2T1 = T2 / mix1.T;
            % Compute STOP criteria
            aux = compute_STOP(x);
            if aux > STOP && it > itMax / 5
                fprintf('NR not converging. It: %d, STOP: %1.2e\n', it, aux);
                return
            end

            STOP = aux;
            % Debug
            % aux_lambda(it) = lambda;
            % aux_STOP(it) = STOP;
        end
        
        % Debug
        % debug_plot_error(it, aux_STOP, aux_lambda);
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

function [p2, T2, p2p1, T2T1] = get_guess(self, mix1, mix2, TN)

    if isempty(mix2)
        M1 = mix1.u / mix1.sound;
        p2p1 = (2 * mix1.gamma * M1^2 - mix1.gamma + 1) / (mix1.gamma + 1);
        T2T1 = p2p1 * (2 / M1^2 + mix1.gamma - 1) / (mix1.gamma + 1);

        if M1 > TN.Mach_thermo
            % Estimate post-shock state considering h2 = h1 + u1^2 / 2
            mix1.h = (enthalpy_mass(mix1) * 1e3 + velocity_relative(mix1)^2/2) * mass(mix1) * 1e-3; % [kJ]
            self.PD.ProblemType = 'HP';
            mix2 = equilibrate(self, mix1, p2p1 * mix1.p);
            T2 = mix2.T; % [K]
        else
            T2 = T2T1 * mix1.T; % [K]
        end

        p2 = p2p1 * convert_bar_to_Pa(mix1.p); % [Pa]
    else
        p2 = convert_bar_to_Pa(mix2.p); % [Pa]
        T2 = mix2.T; % [K]

        p2p1 = p2 / convert_bar_to_Pa(mix1.p);
        T2T1 = T2 / mix1.T;
    end
    
end

function [J, b, guess_moles] = update_system(self, mix1, p2, T2, R0, guess_moles, FLAG_FAST)
    % Update Jacobian matrix and vector b
    r1 = mix1.rho;
    p1 = convert_bar_to_Pa(mix1.p); % [Pa]
    T1 = mix1.T;
    u1 = mix1.u;
    W1 = mix1.W * 1e-3; % [kg/mol]
    h1 = mix1.h / mix1.mi * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [mix2, r2, dVdT_p, dVdp_T] = state(self, mix1, T2, p2, guess_moles);

    W2 = mix2.W * 1e-3;
    h2 = mix2.h / mix2.mi * 1e3; % [J/kg]
    cP2 = mix2.cP / mix2.mi; % [J/(K-kg)]

    alpha = (W1 * u1^2) / (R0 * T1);
    J1 = -r1 / r2 * alpha * dVdp_T - p2 / p1;
    J2 = -r1 / r2 * alpha * dVdT_p;
    b1 = p2 / p1 - 1 + alpha * (r1 / r2 - 1);

    J3 = -u1^2 / R0 * (r1 / r2)^2 * dVdp_T + T2 / W2 * (dVdT_p - 1);
    J4 = -u1^2 / R0 * (r1 / r2)^2 * dVdT_p - T2 * cP2 / R0;
    b2 = (h2 - h1) / R0 - u1^2 / (2 * R0) * (1 - (r1 / r2)^2);

    J = [J1 J2; J3 J4];
    b = [b1; b2];

    % Update guess moles
    if FLAG_FAST
        guess_moles = mix2.Xi * mix2.N;
    end

end

function [mix2, r2, dVdT_p, dVdp_T] = state(self, mix1, T, p, guess_moles)
    % Calculate frozen state given T & p
    self.PD.ProblemType = 'TP';
    p = convert_Pa_to_bar(p); % [bar]
    mix2 = equilibrate_T(self, mix1, p, T, guess_moles);
    r2 = mix2.rho;
    dVdT_p = mix2.dVdT_p;
    dVdp_T = mix2.dVdp_T;
end

function relax = relax_factor(x)
    % Compute relaxation factor
    factor = [0.40546511; 0.40546511];
    lambda = factor ./ abs(x);
    relax = min(1, min(lambda));
end

function [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda)
    % Compute new estimates
    log_p2p1 = log(p2p1) + lambda * x(1);
    log_T2T1 = log(T2T1) + lambda * x(2);
end

function [p2, T2] = apply_antilog(mix1, log_p2p1, log_T2T1)
    % compute p2 and T2
    p2 = exp(log_p2p1) * convert_bar_to_Pa(mix1.p); % [Pa]
    T2 = exp(log_T2T1) * mix1.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function mix2 = save_state(self, mix1, T2, p2, STOP)
    mix2 = state(self, mix1, T2, p2, []);
    mix2.v_shock = mix1.u * mix1.rho / mix2.rho;
    mix2.u = mix1.u - mix2.v_shock; % velocity postshock [m/s] - laboratory fixed
    mix2.error_problem = STOP;
end

function print_convergence(STOP, TOL, T)

    if STOP > TOL
        fprintf('***********************************************************\n')
        fprintf('Convergence error: %1.2e\n', STOP);
    end

    if T > 2e4
        fprintf('***********************************************************\n')
        fprintf('Validity of the results compromise: T = %d K\n', round(T))
        fprintf('Thermodynamic properties fitted to 20000 K\n');
    end

end

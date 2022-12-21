function [mix1, mix2] = det_cj(self, mix1, varargin)
    % Compute pre-shock and post-shock states of a Chapman-Jouguet detonation
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %
    % Optional Args:
    %     mix2 (struct): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing:
    %
    %     * mix1 (struct): Properties of the mixture in the pre-shock state
    %     * mix2 (struct): Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2, guess_moles] = unpack(self, mix1, varargin);
    % Abbreviations
    C = self.C;
    TN = self.TN;
    % Constants
    R0 = C.R0; % Universal gas constant [J/(mol-K)]
    % Definitions
    FLAG_FAST = TN.FLAG_FAST;
    % Solve Chapman-Jouguet detonation
    [T2, p2, STOP, T_guess, p2_guess] = solve_cj_detonation(FLAG_FAST);
    % Check convergence
    print_convergence(STOP, TN.tol_shocks);
    % Save state
    [mix1, mix2] = save_state(self, mix1, T2, p2, STOP);
    mix2.T_guess = T_guess;
    mix2.p_guess = p2_guess;

    % NESTED-FUNCTIONS
    function [T2, p2, STOP, T_guess, p2_guess] = solve_cj_detonation(FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = TN.it_shocks;
        % Initial estimates of p2/p1 and T2/T1
        [p2, T2, p2p1, T2T1, STOP] = get_guess(self, mix1, mix2);
        T_guess = T2; p2_guess = p2;
        % Check FLAG
        if ~FLAG_FAST, guess_moles = []; end
        % Loop
        while STOP > TN.tol_shocks && it < itMax
            % Update iteration
            it = it + 1;
            % Construction of the Jacobian matrix and vector b
            [J, b, guess_moles] = update_system(self, mix1, p2, T2, R0, guess_moles, FLAG_FAST);
            % Solve of the linear system A*x = b
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
            STOP = compute_STOP(x);
        end

    end

end

% SUB-PASS FUNCTIONS
function [self, mix1, mix2, guess_moles] = unpack(self, mix1, x)
    % Unpack input data
    try
        mix2 = x{1};
        guess_moles = mix2.Xi * mix2.N;
    catch
        mix2 = [];
        guess_moles = [];
    end

end

function [p2, T2, p2p1, T2T1, STOP] = get_guess(self, mix1, mix2)

    if isempty(mix2)

        try
            [p2p1, T2T1, ~, ~, Q, STOP] = det_compute_guess(self, mix1, mix1.phi, 1);
            % print_guess(T2T1, p2p1, mix1.T, mix1.p, Q)
        catch
            [p2p1, T2T1, STOP] = det_compute_guess_CEA(self, mix1);
        end

        p2 = p2p1 * convert_bar_to_Pa(mix1.p); % [Pa]
        T2 = T2T1 * mix1.T; % [K]
    else
        p2 = convert_bar_to_Pa(mix2.p); % [Pa]
        T2 = mix2.T; % [K]
        p2p1 = p2 / (convert_bar_to_Pa(mix1.p));
        T2T1 = T2 / mix1.T;
        STOP = 1;
    end

end

function print_guess(T2T1, p2p1, T1, p1, Q)
    fprintf('\n\n START \n\n');
    fprintf('   GUESS:\n');
    fprintf('T2   [K]: %.2f\n', T2T1 * T1);
    fprintf('p2 [bar]: %.2f\n', p2p1 * p1);
    fprintf('Q    [-]: %.2f\n', Q);
    fprintf('\n\n END \n\n');
end

function [J, b, guess_moles] = update_system(self, mix1, p2, T2, R0, guess_moles, FLAG_FAST)
    % Update Jacobian matrix and vector b
    r1 = density(mix1);
    p1 = convert_bar_to_Pa(mix1.p); % [Pa]
    h1 = enthalpy_mass(mix1) * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [mix2, r2, dVdT_p, dVdp_T] = state(self, mix1, T2, p2, guess_moles);

    W2 = MolecularWeight(mix2) * 1e-3; % [kg/mol]
    h2 = enthalpy_mass(mix2) * 1e3; % [J/kg]
    cP2 = cp_mass(mix2) * 1e3; % [J/(K-kg)]
    gamma2_s = adiabaticIndex_sound(mix2); % [-]

    J1 = p1 / p2 + r2 / r1 * gamma2_s * dVdp_T;
    J2 = r2 / r1 * gamma2_s * dVdT_p;
    b1 = p1 / p2 - 1 + gamma2_s * (r2 / r1 - 1);

    J3 = gamma2_s * T2 / (2 * W2) * ((r2 / r1)^2 - 1 - dVdp_T * (1 + (r2 / r1)^2)) ...
        + T2 / W2 * (dVdT_p - 1);
    J4 = -gamma2_s * T2 / (2 * W2) * ((r2 / r1)^2 + 1) * dVdT_p - T2 * cP2 / R0;
    b2 = (h2 - h1) / R0 - gamma2_s * T2 / (2 * W2) * ((r2 / r1)^2 - 1);

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
    r2 = density(mix2);
    dVdT_p = mix2.dVdT_p;
    dVdp_T = mix2.dVdp_T;
end

function relax = relax_factor(x)
    % Compute relaxation factor
    factor = [0.40546511; 0.04879016];
    %     factor = [0.40546511; 0.40546511];
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

function [mix1, mix2] = save_state(self, mix1, T2, p2, STOP)
    mix2 = state(self, mix1, T2, p2, []);
    mix2.u = mix2.sound; % velocity postshock [m/s] - laboratory fixed
    mix2.error_problem = STOP;
    mix1.u = mix2.u * mix2.rho / mix1.rho;
    mix2.v_shock = mix1.u * mix1.rho / mix2.rho;
end

function print_convergence(STOP, TOL)

    if STOP > TOL
        fprintf('Convergence error: %.2f%%\n', STOP);
    end

end

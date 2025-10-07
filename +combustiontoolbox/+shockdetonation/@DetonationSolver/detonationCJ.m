function [mix1, mix2] = detonationCJ(obj, mix1, varargin)
    % Compute pre-shock and post-shock states of a Chapman-Jouguet detonation
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (DetonationSolver): DetonationSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %
    % Optional Args:
    %     mix2 (Mixture): Properties of the mixture in the post-shock state (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     * mix2 (Mixture): Properties of the mixture in the post-shock state
    %
    % Examples:
    %     * [mix1, mix2] = detonationCJ(DetonationSolver(), mix1)
    %     * [mix1, mix2] = detonationCJ(DetonationSolver(), mix1, mix2)

    % Unpack input data
    [mix1, mix2, guess_moles] = unpack(mix1, varargin{:});
    
    % Definitions
    FLAG_FAST = obj.equilibriumSolver.FLAG_FAST;
    R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(mol-K)]
    
    % Solve Chapman-Jouguet detonation
    try
        [T2, p2, STOP, it, T_guess, p2_guess] = solve_cj_detonation(FLAG_FAST);
        assert(STOP < obj.tol0);
    catch
        % If solution has not converged, repeat without composition estimate
        fprintf('Recalculating: %.2f [K]\n', T2);
        guess_moles = [];
        FLAG_FAST = false;
        [T2, p2, STOP, it, T_guess, p2_guess] = solve_cj_detonation(FLAG_FAST);
    end

    % Check convergence
    combustiontoolbox.utils.printConvergence(it, obj.itMax, T2, STOP, obj.tol0);
    
    % Save state
    [mix1, mix2] = save_state(mix1, mix2, STOP);

    % NESTED-FUNCTIONS
    function [T2, p2, STOP, it, T_guess, p2_guess] = solve_cj_detonation(FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = obj.itMax;

        % Initial estimates of p2/p1 and T2/T1
        [p2, T2, p2p1, T2T1, STOP] = get_guess(obj, mix1, mix2);
        T_guess = T2; p2_guess = p2;

        % Check STOP
        if STOP < obj.tol0
            % Set pressure and temperature of mix2
            mix2.p = p2; mix2.T = T2;
        
            % Calculate state given T & p
            mix2 = state(obj.equilibriumSolver, mix1, mix2, []);
        end

        % Check FLAG
        if ~FLAG_FAST, guess_moles = []; end

        % Loop
        while STOP > obj.tol0 && it < itMax
            % Update iteration
            it = it + 1;

            % Construction of the Jacobian matrix and vector b
            [J, b, guess_moles] = update_system(obj.equilibriumSolver, mix1, mix2, p2, T2, R0, guess_moles, FLAG_FAST);

            % Solve of the linear system J*x = b
            x = linsolve(J, b);

            % Calculate correction factor
            lambda = relax_factor(x);

            % Apply correction
            [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda);

            % Apply antilog
            [p2, T2] = apply_antilog(mix1, log_p2p1, log_T2T1); % [bar] and [K]

            % Update ratios
            p2p1 = p2 / mix1.p;
            T2T1 = T2 / mix1.T;
            
            % Compute STOP criteria
            STOP = compute_STOP(x);
        end

    end

end

% SUB-PASS FUNCTIONS
function [mix1, mix2, guess_moles] = unpack(mix1, varargin)
    % Unpack input data

    if nargin > 1
        mix2 = varargin{1}.copy();
        guess_moles = mix2.Xi * mix2.N;
        return
    end

    mix2 = mix1.copy();
    guess_moles = [];
end

function [p2, T2, p2p1, T2T1, STOP] = get_guess(obj, mix1, mix2)

    if mix1.T == mix2.T
        [p2p1, T2T1, STOP] = detonationGuess(obj, mix1);

        p2 = p2p1 * mix1.p; % [bar]
        T2 = T2T1 * mix1.T; % [K]
    else
        p2 = mix2.p; % [bar]
        T2 = mix2.T; % [K]
        p2p1 = p2 / mix1.p;
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

function [J, b, guess_moles] = update_system(equilibriumSolver, mix1, mix2, p2, T2, R0, guess_moles, FLAG_FAST)
    % Update Jacobian matrix and vector b
    r1 = mix1.rho; % [kg/m3]
    p1 = mix1.p; % [bar]
    h1 = mix1.h / mix1.mi; % [J/kg]
    
    % Set pressure and temperature of mix2
    mix2.p = p2; mix2.T = T2;

    % Calculate state given T & p
    [mix2, r2, dVdT_p, dVdp_T] = state(equilibriumSolver, mix1, mix2, guess_moles);

    r2r1 = r2 / r1; % [-]
    p1p2 = p1 / p2; % [-]
    W2 = mix2.W; % [kg/mol]
    h2 = mix2.h / mix2.mi; % [J/kg]
    cp2 = mix2.cp / mix2.mi; % [J/(K-kg)]
    gamma2_s = mix2.gamma_s; % [-]

    J1 = p1p2 + r2r1 * gamma2_s * dVdp_T;
    J2 = r2r1 * gamma2_s * dVdT_p;
    b1 = p1p2 - 1 + gamma2_s * (r2r1 - 1);

    J3 = gamma2_s * T2 / (2 * W2) * (r2r1^2 - 1 - dVdp_T * (1 + r2r1^2)) ...
        + T2 / W2 * (dVdT_p - 1);
    J4 = - gamma2_s * T2 / (2 * W2) * (r2r1^2 + 1) * dVdT_p - T2 * cp2 / R0;
    b2 = (h2 - h1) / R0 - gamma2_s * T2 / (2 * W2) * (r2r1^2 - 1);

    J = [J1 J2; J3 J4];
    b = [b1; b2];

    % Update guess moles
    if FLAG_FAST
        guess_moles = mix2.Xi * mix2.N;
    end

end

function [mix2, r2, dVdT_p, dVdp_T] = state(equilibriumSolver, mix1, mix2, guess_moles)
    % Calculate state given T & p
    equilibriumSolver.problemType = 'TP';
    equilibriumSolver.equilibrateT(mix1, mix2, mix2.T, guess_moles);
    r2 = mix2.rho;
    dVdT_p = mix2.dVdT_p;
    dVdp_T = mix2.dVdp_T;
end

function relax = relax_factor(x)
    % Compute relaxation factor
    factor = 0.40546511;
    lambda = factor ./ abs(x);
    relax = min([1; lambda]);
end

function [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda)
    % Compute new estimates
    log_p2p1 = log(p2p1) + lambda * x(1);
    log_T2T1 = log(T2T1) + lambda * x(2);
end

function [p2, T2] = apply_antilog(mix1, log_p2p1, log_T2T1)
    % compute p2 and T2
    p2 = exp(log_p2p1) * mix1.p; % [bar]
    T2 = exp(log_T2T1) * mix1.T; % [K]
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function [mix1, mix2] = save_state(mix1, mix2, STOP)
    mix2.u = mix2.sound; % velocity postshock [m/s] - laboratory fixed
    mix1.u = mix2.u * mix2.rho / mix1.rho;
    mix1.uShock = mix1.u;
    mix1.cjSpeed = mix1.u;
    mix1.mach = mix1.u / mix1.sound;
    mix2.uShock = mix1.u * mix1.rho / mix2.rho;
    mix2.errorProblem = STOP;
end
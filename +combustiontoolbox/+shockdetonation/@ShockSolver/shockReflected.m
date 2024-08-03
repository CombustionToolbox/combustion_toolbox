function [mix1, mix2, mix5] = shockReflected(obj, mix1, mix2, varargin)
    % Compute pre-shock and post-shock states of a planar reflected shock wave
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     mix2 (Mixture): Properties of the mixture at the post-shock state of the incident shock
    %
    % Optional Args:
    %     mix5 (Mixture): Properties of the mixture in the post-shock state of the reflected shock (previous calculation)
    %
    % Returns:
    %     Tuple containing
    %
    %     * mix1 (Mixture): Properties of the mixture in the pre-shock state of the incident shock
    %     * mix2 (Mixture): Properties of the mixture at the post-shock state of the incident shock
    %     * mix5 (Mixture): Properties of the mixture in the post-shock state of the reflected shock
    %
    % Examples:
    %     * [mix1, mix2, mix5] = shockReflected(ShockSolver(), mix1, mix2)
    %     * [mix1, mix2, mix5] = shockReflected(ShockSolver(), mix1, mix2, mix5)

    % Unpack input data
    [mix5, guess_moles] = unpack(mix2, varargin{:});

    % Definitions
    FLAG_FAST = obj.equilibriumSolver.FLAG_FAST;
    R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(mol-K)]

    % Solve shock reflected
    try
        [T5, p5, STOP, it] = solve_shock_reflected(FLAG_FAST);
        assert(STOP < obj.tol0);
    catch
        % If solution has not converged, repeat without composition estimate
        fprintf('Recalculating: %.2f [m/s]\n', mix1.u);
        guess_moles = [];
        FLAG_FAST = false;
        [T5, p5, STOP, it] = solve_shock_reflected(FLAG_FAST);
    end

    % Check convergence
    combustiontoolbox.utils.printConvergence(it, obj.itMax, T5, STOP, obj.tol0);
    
    % Save state
    mix5 = save_state(mix2, mix5, STOP);

    % NESTED-FUNCTIONS
    function [T5, p5, STOP, it] = solve_shock_reflected(FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = obj.itMax; STOP = 1;
        % Initial estimates of p5/p2 and T5/T2
        [p5, T5, p5p2, T5T2] = get_guess(mix2, mix5);
        % Check FLAG
        if ~FLAG_FAST, guess_moles = []; end
        % Loop
        while STOP > obj.tol0 && it < itMax
            % Update iteration
            it = it + 1;
            % Construction of the Jacobian matrix and vector b
            [J, b, guess_moles] = update_system(obj.equilibriumSolver, mix2, mix5, p5, T5, R0, guess_moles, FLAG_FAST);
            % Solve of the linear system J*x = b
            x = J \ b;
            % Calculate correction factor
            lambda = relax_factor(x);
            % Apply correction
            [log_p5p2, log_T5T2] = apply_correction(x, p5p2, T5T2, lambda);
            % Apply antilog
            [p5, T5] = apply_antilog(mix2, log_p5p2, log_T5T2); % [Pa] and [K]
            % Update ratios
            p5p2 = p5 / (convert_bar_to_Pa(mix2.p));
            T5T2 = T5 / mix2.T;
            % Compute STOP criteria
            STOP = compute_STOP(x);
        end

    end

end

% NESTED FUNCTIONS
function [mix5, guess_moles] = unpack(mix2, varargin)
    % Unpack input data
    try
        mix5 = varargin{1}.copy();
        guess_moles = mix5.Xi * mix5.N;
    catch
        mix5 = mix2.copy();
        guess_moles = [];
    end
    
end

function [p5, T5, p5p2, T5T2] = get_guess(mix2, mix5)

    if mix2.u == mix5.u
        T5T2 = 2;

        b = (mix2.gamma_s + 1) / (mix2.gamma_s - 1);
        p5p2 = 0.5 * (b + sqrt(8 + b^2)); % positive root of (7.45 - report CEA)

        p5 = p5p2 * convert_bar_to_Pa(mix2.p); % [Pa]
        T5 = T5T2 * mix2.T; % [K]
    else
        p5 = convert_bar_to_Pa(mix5.p); % [Pa]
        T5 = mix5.T; % [K]

        p5p2 = p5 / (convert_bar_to_Pa(mix2.p));
        T5T2 = T5 / mix2.T;
    end

end

function [J, b, guess_moles] = update_system(equilibriumSolver, mix2, mix5, p5, T5, R0, guess_moles, FLAG_FAST)
    % Update Jacobian matrix and vector b
    r2 = mix2.rho; % [kg/m3]
    p2 = convert_bar_to_Pa(mix2.p); % [Pa]
    T2 = mix2.T; % [K]
    u2 = mix2.u; % [m/s]
    W2 = mix2.W; % [kg/mol]
    h2 = mix2.h / mix2.mi; % [J/kg]

    % Set pressure and temperature of mix5
    mix5.p = convert_Pa_to_bar(p5); mix5.T = T5;

    % Calculate state given T & p
    [mix5, r5, dVdT_p, dVdp_T] = state(equilibriumSolver, mix2, mix5, guess_moles);

    W5 = mix5.W; % [kg/mol]
    h5 = mix5.h / mix5.mi; % [J/kg]
    cp5 = mix5.cp / mix5.mi; % [J/(K-kg)]

    alpha = (W2 * u2^2) / (R0 * T2);
    J1 = (r5 / r2) / (r5 / r2 - 1)^2 * alpha * dVdp_T - p5 / p2;
    J2 = (r5 / r2) / (r5 / r2 - 1)^2 * alpha * dVdT_p;
    b1 = p5 / p2 - 1 - alpha * r5 / r2 / (r5 / r2 - 1);

    J3 = u2^2 / R0 * (r5 / r2) / (r5 / r2 - 1)^2 * dVdp_T + T5 / W5 * (dVdT_p - 1);
    J4 = u2^2 / R0 * (r5 / r2) / (r5 / r2 - 1)^2 * dVdT_p - T5 * cp5 / R0;
    b2 = (h5 - h2) / R0 - u2^2 / (2 * R0) * (r5 / r2 + 1) / (r5 / r2 - 1);

    J = [J1 J2; J3 J4];
    b = [b1; b2];

    % Update guess moles
    if FLAG_FAST
        guess_moles = mix5.Xi * mix5.N;
    end

end

function [mix5, r5, dVdT_p, dVdp_T] = state(equilibriumSolver, mix2, mix5, guess_moles)
    % Calculate state given T & p
    equilibriumSolver.problemType = 'TP';
    equilibrateT(equilibriumSolver, mix2, mix5, mix5.T, guess_moles);
    r5 = mix5.rho;
    dVdT_p = mix5.dVdT_p;
    dVdp_T = mix5.dVdp_T;
end

function relax = relax_factor(x)
    % Compute relaxation factor
    factor = [0.40546511; 0.40546511];
    lambda = factor ./ abs(x);
    relax = min(1, min(lambda));
end

function [log_p5p2, log_T5T2] = apply_correction(x, p5p2, T5T2, lambda)
    % Compute new estimates
    log_p5p2 = log(p5p2) + lambda * x(1);
    log_T5T2 = log(T5T2) + lambda * x(2);
end

function [p5, T5] = apply_antilog(mix2, log_p5p2, log_T5T2)
    % compute p2 and T2
    p5 = exp(log_p5p2) * convert_bar_to_Pa(mix2.p); % [Pa]
    T5 = exp(log_T5T2) * mix2.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function mix5 = save_state(mix2, mix5, STOP)
    mix5.u = convert_bar_to_Pa(mix5.p - mix2.p) / (mix2.u * mix2.rho) - mix2.u;
    mix5.uShock = mix2.u * mix2.rho / mix5.rho;
    mix5.mach = mix5.uShock / mix5.sound;
    mix5.errorProblem = STOP;
end
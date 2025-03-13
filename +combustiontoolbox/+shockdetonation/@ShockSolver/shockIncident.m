function [mix1, mix2] = shockIncident(obj, mix1, u1, varargin)
    % Compute pre-shock and post-shock states of a planar incident shock wave
    %
    % This method is based on the method outlined in Gordon, S., & McBride,
    % B. J. (1994). NASA reference publication, 1311.
    %
    % Args:
    %     obj (ShockSolver): ShockSolver object
    %     mix1 (Mixture): Properties of the mixture in the pre-shock state
    %     u1 (float): Pre-shock velocity [m/s]
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
    %     * [mix1, mix2] = shockIncident(ShockSolver(), mix1, u1)
    %     * [mix1, mix2] = shockIncident(ShockSolver(), mix1, u1, mix2)

    % Unpack input data
    [mix1, mix2, guess_moles] = unpack(mix1, u1, varargin{:});

    if obj.equilibriumSolver.FLAG_TCHEM_FROZEN
        STOP = 0;
        [~, p2p1, T2T1, ~, ~] = obj.shockIncidentIdeal(mix1.gamma, mix1.mach);
        T2 = T2T1 * mix1.T; % [K]
        p2 = p2p1 * mix1.p; % [bar]
        mix2.p = p2; mix2.T = T2;
        mix2 = equilibrateT(obj.equilibriumSolver, mix1, mix2, T2);
        mix2 = save_state(mix1, mix2, STOP);
        return
    end

    % Definitions
    FLAG_FAST = obj.equilibriumSolver.FLAG_FAST;
    R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(mol-K)]

    % Solve shock incident
    try
        [T2, p2, STOP, it] = solve_shock_incident(FLAG_FAST); % p2 [Pa]
        assert(STOP < obj.tol0);
    catch
        % If solution has not converged, repeat without composition estimate
        fprintf('Recalculating: %.2f [m/s]\n', u1);
        guess_moles = [];
        FLAG_FAST = false;
        [T2, p2, STOP, it] = solve_shock_incident(FLAG_FAST); % p2 [Pa]
    end

    % Check convergence
    combustiontoolbox.utils.printConvergence(it, obj.itMax, T2, STOP, obj.tol0);
    
    % Save state
    mix2 = save_state(mix1, mix2, STOP);

    % NESTED-FUNCTIONS
    function [T2, p2, STOP, it] = solve_shock_incident(FLAG_FAST)
        % Miscellaneous
        it = 0; itMax = obj.itMax; STOP = 1;
        % Initial estimates of p2/p1 and T2/T1
        [p2, T2, p2p1, T2T1] = get_guess(obj, mix1, mix2);
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
function [mix1, mix2, guess_moles] = unpack(mix1, u1, varargin)
    % Unpack input data
    mix1.u = u1; % pre-shock velocity [m/s] - laboratory fixed
    mix1.uShock = u1; % pre-shock velocity [m/s] - shock fixed
    mix1.mach = mix1.u / mix1.sound; % pre-shock Mach number [-]

    try
        mix2 = varargin{1}.copy();
        guess_moles = mix2.Xi * mix2.N;
    catch
        mix2 = mix1.copy();
        guess_moles = [];
    end
    
end

function [p2, T2, p2p1, T2T1] = get_guess(obj, mix1, mix2)

    if mix1.u == mix2.u
        M1 = mix1.u / mix1.sound;
        p2p1 = (2 * mix1.gamma * M1^2 - mix1.gamma + 1) / (mix1.gamma + 1);
        T2T1 = p2p1 * (2 / M1^2 + mix1.gamma - 1) / (mix1.gamma + 1);

        if M1 > obj.machThermo
            % Estimate post-shock state considering h2 = h1 + u1^2 / 2
            mix2.h = mix1.h + 0.5 * mix1.u^2 * mix1.mi; % [J]

            % Initialize mix2 and set pressure
            mix2.p = p2p1 * mix1.p;
            
            % Solve chemical transformation
            obj.equilibriumSolver.problemType = 'HP';
            solve(obj.equilibriumSolver, mix2);

            % Get temperature
            T2 = mix2.T; % [K]
        else
            T2 = T2T1 * mix1.T; % [K]
        end

        p2 = p2p1 * mix1.p; % [bar]
    else
        p2 = mix2.p; % [bar]
        T2 = mix2.T; % [K]

        p2p1 = p2 / mix1.p;
        T2T1 = T2 / mix1.T;
    end
    
end

function [J, b, guess_moles] = update_system(equilibriumSolver, mix1, mix2, p2, T2, R0, guess_moles, FLAG_FAST)
    % Update Jacobian matrix and vector b
    r1 = mix1.rho; % [kg/m3]
    p1 = mix1.p; % [bar]
    T1 = mix1.T; % [K]
    u1 = mix1.u; % [m/s]
    W1 = mix1.W; % [kg/mol]
    h1 = mix1.h / mix1.mi; % [J/kg]
    
    % Set pressure and temperature of mix2
    mix2.p = p2; mix2.T = T2;

    % Calculate state given T & p
    [mix2, r2, dVdT_p, dVdp_T] = state(equilibriumSolver, mix1, mix2, guess_moles);
    
    r1r2 = r1 / r2; % [-]
    p2p1 = p2 / p1; % [-]
    W2 = mix2.W; % [kg/mol]
    h2 = mix2.h / mix2.mi; % [J/kg]
    cp2 = mix2.cp / mix2.mi; % [J/(K-kg)]

    alpha = (W1 * u1^2) / (R0 * T1);
    J1 = -r1r2 * alpha * dVdp_T - p2p1;
    J2 = -r1r2 * alpha * dVdT_p;
    b1 = p2p1 - 1 + alpha * (r1 / r2 - 1);

    J3 = -u1^2 / R0 * r1r2^2 * dVdp_T + T2 / W2 * (dVdT_p - 1);
    J4 = -u1^2 / R0 * r1r2^2 * dVdT_p - T2 * cp2 / R0;
    b2 = (h2 - h1) / R0 - u1^2 / (2 * R0) * (1 - r1r2^2);

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

function mix2 = save_state(mix1, mix2, STOP)
    mix2.uShock = mix1.u * mix1.rho / mix2.rho;
    mix2.u = mix1.u - mix2.uShock; % velocity postshock [m/s] - laboratory fixed
    mix2.mach = mix2.uShock / mix2.sound;
    mix2.errorProblem = STOP;
end
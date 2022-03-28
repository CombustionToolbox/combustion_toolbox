function [mix1, mix2, mix5] = shock_reflected(self, mix1, u1, mix2, varargin)
    % Compute pre-shock and post-shock states of a planar reflected shock wave
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the mixture in the pre-shock state of the incident shock
    %     u1 (float):    Pre-shock velocity [m/s]
    %     mix2 (struct): Properties of the mixture at the post-shock state of the incident shock
    %
    % Optional Args:
    %     mix5 (struct): Properties of the mixture in the post-shock state of the reflected shock (previous calculation)
    %
    % Returns:
    %     mix1 (struct): Properties of the mixture in the pre-shock state of the incident shock
    %     mix2 (struct): Properties of the mixture at the post-shock state of the incident shock
    %     mix5 (struct): Properties of the mixture in the post-shock state of the reflected shock

    % Unpack input data
    [self, mix1, mix2, mix5] = unpack(self, mix1, u1, mix2, varargin);
    % Abbreviations
    C = self.C;
    TN = self.TN;
    % Constants
    R0 = C.R0; % Universal gas constant [J/(mol-K)]
    % Miscellaneous
    it = 0;
    itMax = TN.it_shocks;
    STOP = 1.;
    % Initial estimates of p5/p2 and T5/T2
    [p5, T5, p5p2, T5T2] = get_guess(mix2, mix5);
    % Loop
    while STOP > TN.tol_shocks && it < itMax
        it = it + 1;
        % Construction of the Jacobian matrix and vector b
        [J, b] = update_system(self, mix2, p5, T5, R0);
        % Solve of the linear system A*x = b
        x = J\b;
        % Calculate correction factor
        lambda = relax_factor(x);
        % Apply correction
        [log_p5p2, log_T5T2] = apply_correction(x, p5p2, T5T2, lambda);
        % Apply antilog
        [p5, T5] = apply_antilog(mix2, log_p5p2, log_T5T2); % [Pa] and [K]
        % Update ratios
        p5p2 = p5 / (mix2.p * 1e5);
        T5T2 = T5 / mix2.T;
        % Compute STOP criteria
        STOP = compute_STOP(x);
    end
    % Check convergence
    print_convergence(STOP, TN.tol_shocks);
    % Save state
    mix5 = save_state(self, mix2, T5, p5, STOP);
end
% NESTED FUNCTIONS
function [self, mix1, mix2, mix5] = unpack(self, mix1, u1, mix2, x)
    % Unpack input data
    if length(x) > 0
        mix5 = x{1};
    else
        mix5 = [];
    end
end

function [p5, T5, p5p2, T5T2] = get_guess(mix2, mix5)
    if isempty(mix5)
        T5T2 = 2;
        aux  = roots([1, -(mix2.gamma_s + 1) / (mix2.gamma_s - 1), -2]); % roots of (7.45 - report CEA)
        p5p2 = aux(aux > 0); % positive root

        p5 = p5p2 * mix2.p * 1e5; % [Pa]
        T5 = T5T2 * mix2.T;       % [K]
    else
        p5 = mix5.p * 1e5; % [Pa]
        T5 = mix5.T;       % [K]

        p5p2 = p5 / (mix2.p * 1e5);
        T5T2 = T5 / mix2.T;
    end
end

function [J, b] = update_system(self, mix2, p5, T5, R0)
    % Update Jacobian matrix and vector b
    r2 = mix2.rho;
    p2 = mix2.p *1e5; % [Pa]
    T2 = mix2.T;
    u2 = mix2.u;
    W2 = mix2.W * 1e-3; % [kg/mol]
    h2 = mix2.h / mix2.mi * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [mix5, r5, dVdT_p, dVdp_T] = state(self, mix2, T5, p5);
    
    W5 = mix5.W * 1e-3;
    h5 = mix5.h / mix5.mi * 1e3; % [J/kg]
    cP5 = mix5.cP / mix5.mi; % [J/(K-kg)]
    
    alpha = (W2 * u2^2) / (R0 * T2);
    J1 = (r5/r2) / (r5/r2 - 1)^2 * alpha * dVdp_T - p5 / p2;
    J2 = (r5/r2) / (r5/r2 - 1)^2 * alpha * dVdT_p;
    b1 = p5/p2 - 1 - alpha * r5/r2 / (r5/r2 - 1);
    
    J3 = u2^2 / R0 * (r5/r2) / (r5/r2 - 1)^2 * dVdp_T + T5 / W5 * (dVdT_p - 1);
    J4 = u2^2 / R0 * (r5/r2) / (r5/r2 - 1)^2 * dVdT_p - T5 * cP5 / R0;
    b2 = (h5 - h2) / R0 - u2^2 / (2*R0) * (r5/r2 + 1) / (r5/r2 - 1);
    
    J = [J1 J2; J3 J4];
    b = [b1; b2];
end

function [mix5, r5, dVdT_p, dVdp_T]= state(self, mix2, T, p)
    % Calculate frozen state given T & p
    self.PD.ProblemType = 'TP';
    p = p*1e-5; % [bar]
    mix5 = equilibrate_T(self, mix2, p, T);
    r5 = mix5.rho;
    dVdT_p = mix5.dVdT_p;
    dVdp_T = mix5.dVdp_T;
end

function relax = relax_factor(x)
    % Compute relaxation factor
    factor = [0.40546511; 0.04879016];
%     factor = [0.40546511; 0.40546511];
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
    p5 = exp(log_p5p2) * mix2.p * 1e5; % [Pa]
    T5 = exp(log_T5T2) * mix2.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function mix5 = save_state(self, mix2, T5, p5, STOP)
    mix5 = state(self, mix2, T5, p5);
    mix5.u = (mix5.p - mix2.p)*1e5 / (mix2.u * mix2.rho) - mix2.u;
    mix5.error_problem = STOP;
    mix5.v_shock = mix2.u * mix2.rho / mix5.rho;
end

function print_convergence(STOP, TOL)
    if STOP > TOL
        fprintf('Convergence error: %.2f\n', STOP);
    end
end


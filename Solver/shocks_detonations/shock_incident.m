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
    %     mix1 (struct): Properties of the mixture in the pre-shock state
    %     mix2 (struct): Properties of the mixture in the post-shock state

    % Unpack input data
    [self, mix1, mix2] = unpack(self, mix1, u1, varargin);
    % Abbreviations 
    C = self.C;
    TN = self.TN;
    % Constants
    R0 = C.R0; % Universal gas constant [J/(mol-K)]
    % Miscellaneous
    it = 0;
    itMax = TN.it_shocks;
    STOP = 1.;
    % Initial estimates of p2/p1 and T2/T1
    [p2, T2, p2p1, T2T1] = get_guess(mix1, mix2, TN);
    % Loop
    while STOP > TN.tol_shocks && it < itMax
        it = it + 1;
        % Construction of the Jacobian matrix and vector b
        [J, b] = update_system(self, mix1, p2, T2, R0);
        % Solve of the linear system A*x = b
        x = J\b;
        % Calculate correction factor
        lambda = relax_factor(x);
        % Apply correction
        [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda);
        % Apply antilog
        [p2, T2] = apply_antilog(mix1, log_p2p1, log_T2T1); % [Pa] and [K]
        % Update ratios
        p2p1 = p2 / (mix1.p * 1e5);
        T2T1 = T2 / mix1.T;
        % Compute STOP criteria
        STOP = compute_STOP(x);
    end
    % Check convergence
    print_convergence(STOP, TN.tol_shocks, T2);
    % Save state
    mix2 = save_state(self, mix1, T2, p2, STOP);
end
% SUB-PASS FUNCTIONS
function [self, mix1, mix2] = unpack(self, mix1, u1, x)
    % Unpack input data
    mix1.u = u1;       % velocity pre-shock [m/s] - laboratory fixed
    mix1.v_shock = u1; % velocity pre-shock [m/s] - shock fixed
    if length(x) > 0
        mix2 = x{1};
    else
        mix2 = [];
    end
end

function [p2, T2, p2p1, T2T1] = get_guess(mix1, mix2, TN)
    if isempty(mix2)
        V1 = 1/mix1.rho; % [m3/kg]
        V = V1/TN.volumeBoundRation;

        % p2p1 = (2*gamma1 * M1^2 - gamma1 + 1) / (gamma1 + 1);
        % T2T1 = p2p1 * (2/M1^2 + gamma1 - 1) / (gamma1 + 1);

        p2p1 = 1 + (mix1.rho * mix1.u^2  / mix1.p * (1 - V/V1)) * 1e-5;
        T2T1 = p2p1 * V / V1;

        p2 = p2p1 * mix1.p * 1e5; % [Pa]
        T2 = T2T1 * mix1.T;       % [K]
           
%         p2p1 = p2/mix1.p;
%         T2T1 = T2/mix1.T;
    else
        p2 = mix2.p * 1e5; % [Pa]
        T2 = mix2.T;       % [K]

        p2p1 = p2 / (mix1.p * 1e5);
        T2T1 = T2 / mix1.T;
    end
end

function [J, b] = update_system(self, mix1, p2, T2, R0)
    % Update Jacobian matrix and vector b
    r1 = mix1.rho;
    p1 = mix1.p *1e5; % [Pa]
    T1 = mix1.T;
    u1 = mix1.u;
    W1 = mix1.W * 1e-3; % [kg/mol]
    h1 = mix1.h / mix1.mi * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [mix2, r2, dVdT_p, dVdp_T] = state(self, mix1, T2, p2);
    
    W2 = mix2.W * 1e-3;
    h2 = mix2.h / mix2.mi * 1e3; % [J/kg]
    cP2 = mix2.cP / mix2.mi; % [J/(K-kg)]
    
    alpha = (W1 * u1^2) / (R0 * T1);
    J1 = -r1/r2 * alpha * dVdp_T - p2 / p1;
    J2 = -r1/r2 * alpha * dVdT_p;
    b1 = p2/p1 - 1 + alpha * (r1/r2 - 1);
    
    J3 = -u1^2 / R0 * (r1/r2)^2 * dVdp_T + T2 / W2 * (dVdT_p - 1);
    J4 = -u1^2 / R0 * (r1/r2)^2 * dVdT_p - T2 * cP2 / R0;
    b2 = (h2 - h1) / R0 - u1^2 / (2*R0) * (1 - (r1/r2)^2);
    
    J = [J1 J2; J3 J4];
    b = [b1; b2];
end

function [mix2, r2, dVdT_p, dVdp_T]= state(self, mix1, T, p)
    % Calculate frozen state given T & p
    self.PD.ProblemType = 'TP';
    p = p*1e-5; % [bar]
    mix2 = equilibrate_T(self, mix1, p, T);
    r2 = mix2.rho;
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
    p2 = exp(log_p2p1) * mix1.p * 1e5; % [Pa]
    T2 = exp(log_T2T1) * mix1.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function mix2 = save_state(self, mix1, T2, p2, STOP)
    mix2 = state(self, mix1, T2, p2);
    mix2.v_shock = mix1.u * mix1.rho / mix2.rho;
    mix2.u = mix1.u - mix2.v_shock; % velocity postshock [m/s] - laboratory fixed
    mix2.error_problem = STOP;
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

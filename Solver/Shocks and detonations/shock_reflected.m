function [str1, str2, str5] = shock_reflected(varargin)
    % Solve planar reflected shock wave

    % Unpack input data
    [self, str1, str2, str5] = unpack(varargin);
    % Abbreviations
    C = self.C;
    TN = self.TN;
    % Constants
    R0 = C.R0; % Universal gas constant [J/(mol-K)]
    % Miscelaneous
    it = 0;
    itMax = TN.it_shocks;
    STOP = 1.;
    % Initial estimates of p5/p2 and T5/T2
    [p5, T5, p5p2, T5T2] = get_guess(str2, str5);
    % Loop
    while STOP > TN.tol_shocks && it < itMax
        it = it + 1;
        % Construction of the Jacobian matrix and vector b
        [J, b] = update_system(self, str2, p5, T5, R0);
        % Solve of the linear system A*x = b
        x = J\b;
        % Calculate correction factor
        lambda = relax_factor(x);
        % Apply correction
        [log_p5p2, log_T5T2] = apply_correction(x, p5p2, T5T2, lambda);
        % Apply antilog
        [p5, T5] = apply_antilog(str2, log_p5p2, log_T5T2); % [Pa] and [K]
        % Update ratios
        p5p2 = p5 / (str2.p * 1e5);
        T5T2 = T5 / str2.T;
        % Compute STOP criteria
        STOP = compute_STOP(x);
    end
    % Check convergence
    print_convergence(STOP, TN.tol_shocks);
    % Save state
    str5 = save_state(self, str2, T5, p5, STOP);
end
% NESTED FUNCTIONS
function [self, str1, str2, str5] = unpack(x)
    % Unpack input data
    self = x{1};
    str1 = x{2};
    u1   = x{3};
    str1.u = u1; % velocity preshock [m/s]
    str2 = x{4};
    if length(x) > 4
        str5 = x{5};
    else
        str5 = [];
    end
end

function [p5, T5, p5p2, T5T2] = get_guess(str2, str5)
    if isempty(str5)
        T5T2 = 2;
        aux  = roots([1, -(str2.gamma_s + 1) / (str2.gamma_s - 1), -2]); % roots of (7.45 - report CEA)
        p5p2 = aux(aux > 0); % positive root

        p5 = p5p2 * str2.p * 1e5; % [Pa]
        T5 = T5T2 * str2.T;       % [K]
    else
        p5 = str5.p * 1e5; % [Pa]
        T5 = str5.T;       % [K]

        p5p2 = p5 / (str2.p * 1e5);
        T5T2 = T5 / str2.T;
    end
end

function [J, b] = update_system(self, str2, p5, T5, R0)
    % Update Jacobian matrix and vector b
    r2 = str2.rho;
    p2 = str2.p *1e5; % [Pa]
    T2 = str2.T;
    u2 = str2.u;
    W2 = str2.W * 1e-3; % [kg/mol]
    h2 = str2.h / str2.mi * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [str5, r5, dVdT_p, dVdp_T] = state(self, str2, T5, p5);
    
    W5 = str5.W * 1e-3;
    h5 = str5.h / str5.mi * 1e3; % [J/kg]
    cP5 = str5.cP / str5.mi; % [J/(K-kg)]
    
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

function [str5, r5, dVdT_p, dVdp_T]= state(self, str2, T, p)
    % Calculate frozen state given T & p
    self.PD.ProblemType = 'TP';
    p = p*1e-5; % [bar]
    str5 = equilibrate_T(self, str2, p, T);
    r5 = str5.rho;
    dVdT_p = str5.dVdT_p;
    dVdp_T = str5.dVdp_T;
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

function [p5, T5] = apply_antilog(str2, log_p5p2, log_T5T2)
    % compute p2 and T2
    p5 = exp(log_p5p2) * str2.p * 1e5; % [Pa]
    T5 = exp(log_T5T2) * str2.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function str5 = save_state(self, str2, T5, p5, STOP)
    str5 = state(self, str2, T5, p5);
    str5.u = (str5.p - str2.p)*1e5 / (str2.u * str2.rho) - str2.u;
    str5.error_problem = STOP;
    str5.v_shock = str2.u * str2.rho / str5.rho;
end

function print_convergence(STOP, TOL)
    if STOP > TOL
        fprintf('Convergence error: %.2f\n', STOP);
    end
end


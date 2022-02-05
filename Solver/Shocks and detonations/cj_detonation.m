function [str1, str2] = cj_detonation(varargin)
% Unpack input data
[self, str1, str2] = unpack(varargin);
% Abbreviations 
C = self.C;
TN = self.TN;
% Constants
R0 = C.R0; % Universal gas constant [J/(mol-K)]
% Miscelaneous
it = 0;
itMax = TN.it_shocks;
% Initial estimates of p2/p1 and T2/T1
[p2, T2, p2p1, T2T1, STOP] = get_guess(self, str1, str2);
T_guess = T2;
p_guess = p2;
% Loop
while STOP > TN.tol_shocks && it < itMax
    it = it + 1;
    % Construction of the Jacobian matrix and vector b
    [J, b] = update_system(self, str1, p2, T2, R0);
    % Solve of the linear system A*x = b
    x = J\b;
    % Calculate correction factor
    lambda = relax_factor(x);
    % Apply correction
    [log_p2p1, log_T2T1] = apply_correction(x, p2p1, T2T1, lambda);
    % Apply antilog
    [p2, T2] = apply_antilog(str1, log_p2p1, log_T2T1); % [Pa] and [K]
    % Update ratios
    p2p1 = p2 / (str1.p * 1e5);
    T2T1 = T2 / str1.T;
    % Compute STOP criteria
    STOP = compute_STOP(x);
end
% Check convergence
print_convergence(STOP, TN.tol_shocks);
% Save state
[str1, str2] = save_state(self, str1, T2, p2, STOP);
str2.T_guess = T_guess;
str2.p_guess = p_guess;
end
% SUB-PASS FUNCTIONS
function [self, str1, str2] = unpack(x)
    % Unpack input data
    self = x{1};
    str1 = x{2};
    if length(x) > 2
        str2 = x{3};
    else
        str2 = [];
    end
end

function [p2, T2, p2p1, T2T1, STOP] = get_guess(self, str1, str2)
    if isempty(str2)
        [p2p1, T2T1, ~, ~, Q, STOP] = compute_guess_det(self, str1, str1.phi, 1);
%             print_guess(T2T1, p2p1, str1.T, str1.p, Q)

        p2 = p2p1 * str1.p * 1e5; % [Pa]
        T2 = T2T1 * str1.T;       % [K]
    else
        p2 = str2.p * 1e5; % [Pa]
        T2 = str2.T;       % [K]
        p2p1 = p2 / (str1.p * 1e5);
        T2T1 = T2 / str1.T;
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

function [J, b] = update_system(self, str1, p2, T2, R0)
    % Update Jacobian matrix and vector b
    r1 = str1.rho;
    p1 = str1.p *1e5; % [Pa]
    h1 = str1.h / str1.mi * 1e3; % [J/kg]
    % Calculate frozen state given T & p
    [str2, r2, dVdT_p, dVdp_T] = state(self, str1, T2, p2);
    
    W2 = str2.W * 1e-3;
    h2 = str2.h / str2.mi * 1e3; % [J/kg]
    cP2 = str2.cP / str2.mi; % [J/(K-kg)]
    gamma2_s = str2.gamma_s;
    
    J1 =  p1/p2 + r2/r1 * gamma2_s * dVdp_T;
    J2 =          r2/r1 * gamma2_s * dVdT_p;
    b1 = p1/p2 - 1 + gamma2_s * (r2/r1 - 1);
    
    J3 =  gamma2_s * T2 / (2*W2) * ((r2/r1)^2 - 1 - dVdp_T * (1 + (r2/r1)^2)) ...
         + T2 / W2 * (dVdT_p - 1);
    J4 = -gamma2_s * T2 / (2*W2) * ((r2/r1)^2 + 1) * dVdT_p - T2 * cP2 / R0;
    b2 = (h2 - h1) / R0 - gamma2_s * T2 / (2*W2) * ((r2/r1)^2 - 1);
    
    J = [J1 J2; J3 J4];
    b = [b1; b2];
end

function [str2, r2, dVdT_p, dVdp_T]= state(self, str1, T, p)
    % Calculate frozen state given T & p
    self.PD.ProblemType = 'TP';
    p = p*1e-5; % [bar]
    str2 = equilibrate_T(self, str1, p, T);
    r2 = str2.rho;
    dVdT_p = str2.dVdT_p;
    dVdp_T = str2.dVdp_T;
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

function [p2, T2] = apply_antilog(str1, log_p2p1, log_T2T1)
    % compute p2 and T2
    p2 = exp(log_p2p1) * str1.p * 1e5; % [Pa]
    T2 = exp(log_T2T1) * str1.T;
end

function STOP = compute_STOP(x)
    % compute stop condition
    STOP = max(abs(x));
end

function [str1, str2] = save_state(self, str1, T2, p2, STOP)
    str2 = state(self, str1, T2, p2);
    str2.u = str2.sound; % velocity postshock [m/s] - laboratory fixed
    str2.error_problem = STOP;
    str1.u = str2.u * str2.rho / str1.rho;
    str2.v_shock = str1.u * str1.rho / str2.rho;
end

function print_convergence(STOP, TOL)
    if STOP > TOL
        fprintf('Convergence error: %.2f%%\n', STOP);
    end
end

function [N0, STOP] = equilibrium_reduced(self, pP, TP, mix1, guess_moles)
    % Obtain equilibrium composition [moles] for the given temperature [K] and pressure [bar].
    % The code stems from the minimization of the free energy of the system by using Lagrange
    % multipliers combined with a Newton-Raphson method, upon condition that initial gas
    % properties are defined by two functions of states. e.g., temperature and pressure.
    % The algorithm implemented take advantage of the sparseness of the
    % upper left submatrix obtaining a matrix A of size NE + NS-NG + 1. 
    %
    % This method is based on Gordon, S., & McBride, B. J. (1994). NASA reference publication,
    % 1311.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     pP (float): Pressure [bar]
    %     TP (float): Temperature [K]
    %     mix1 (struct): Properties of the initial mixture
    %     guess_moles (float): mixture composition [mol] of a previous computation
    %
    % Returns:
    %     Tuple containing
    %
    %     - N0 (float): Equilibrium composition [moles] for the given temperature [K] and pressure [bar]
    %     - STOP (float): Relative error [-] 

    % Generalized Gibbs minimization method
    
    % Abbreviations ---------------------
    S = self.S;
    C = self.C;
    TN = self.TN;
    % -----------------------------------
    
    N0 = C.N0.value;
    A0 = C.A0.value;
    R0TP = C.R0 * TP; % [J/(mol)]
    % Initialization
    NatomE = mix1.NatomE;
    NatomE_tol = NatomE(NatomE > TN.tolE);
    max_NatomE = max(NatomE);
    NP_0 = 0.1;
    NP = NP_0;
    SIZE = -log(TN.tolN);
    
    % Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
    % for the sum of elements whose value is <= tolN
    ind_A0_E0 = remove_elements(NatomE, A0, TN.tolN);
    % List of indices with nonzero values
    [temp_ind_nswt, temp_ind_swt, temp_ind_cryogenic, temp_ind_E, temp_NE] = temp_values(S, NatomE, TN.tolN);
    % Update temp values
    [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, ~] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, NP, SIZE);
    [temp_ind, temp_ind_swt, temp_NS] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic);
    
    temp_ind_swt_0 = temp_ind_swt;
    temp_ind_swt = [];
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NS = length(temp_ind);
    % Initialize species vector N0
    if isempty(guess_moles)
        N0(temp_ind_nswt, 1) = NP_0/temp_NG;
    else
        N0(temp_ind_nswt, 1) = guess_moles(temp_ind_nswt);
        NP_0 = sum(guess_moles);
        NP = NP_0;
    end
    % Standard Gibbs free energy
    g0 = set_g0(S.LS, TP, self.DB);
    % Dimensionless Chemical potential
    muRT = g0/R0TP;
    % Definitions
    A0_T = A0';
    % Solve system
    x = equilibrium_loop;
    % Check condensed species
    [temp_ind, temp_ind_swt, FLAG] = check_condensed_species(A0, x, temp_ind_nswt, temp_ind_swt_0, temp_ind_E, temp_NS, muRT);
    if FLAG
        if any(isnan(x))
            NP = NP_0;
            [temp_ind_nswt, temp_ind_swt, temp_ind_cryogenic, temp_ind_E, temp_NE] = temp_values(S, NatomE, TN.tolN);
            % Update temp values
            [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, ~] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, NP, SIZE);
            [temp_ind, temp_ind_swt, temp_NS] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic);
            % Initialize species vector N0 
            N0(:, 1) = NP_0/temp_NS;
        end
        STOP = 1;
        temp_NS = length(temp_ind);
        equilibrium_loop;
    end
    % NESTED FUNCTION
    function x = equilibrium_loop
        it = 0;
        itMax = 50 + round(S.NS/2);
        STOP = 1.;
        while STOP > TN.tolN && it < itMax
            it = it + 1;
            % Chemical potential
            muRT(temp_ind_nswt) =  g0(temp_ind_nswt) / R0TP + log(N0(temp_ind_nswt, 1) / NP) + log(pP);
            % Construction of matrix A
            A = update_matrix_A(A0_T, N0, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS, temp_NE);
            % Construction of vector b      
            b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, muRT);
            % Solve of the linear system A*x = b
            x = A\b;
            % Extract solution
            pi_i = x(1:temp_NE);
            Delta_nj = x(temp_NE+1:end-1);
            Delta_NP = x(end);
            % Compute correction moles of gases
            Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, Delta_NP, muRT, temp_ind_nswt, temp_ind_E);
            % Calculate correction factor
            lambda = relax_factor(NP, N0(temp_ind, 1), [Delta_ln_nj; Delta_nj], Delta_NP, temp_NG, SIZE);
            % Apply correction
            N0(temp_ind_nswt, 1) = log(N0(temp_ind_nswt, 1)) + lambda * Delta_ln_nj;
            N0(temp_ind_swt, 1) = N0(temp_ind_swt, 1) + lambda * Delta_nj;
            NP_log = log(NP) + lambda * Delta_NP;
            % Apply antilog
            [N0, NP] = apply_antilog(N0, NP_log, temp_ind_nswt);
            % Compute STOP criteria
            STOP = compute_STOP(NP_0, NP, Delta_NP, N0(temp_ind, 1), [Delta_ln_nj; Delta_nj], temp_NG, A0(temp_ind, temp_ind_E), NatomE_tol, max_NatomE, TN.tolE);
            % Update temp values in order to remove species with moles < tolerance
            [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, NP, SIZE);
        end
    end
end

% SUB-PASS FUNCTIONS
function [temp_ind, temp_ind_swt, FLAG] = check_condensed_species(A0, x, temp_ind_nswt, temp_ind_swt_0, temp_ind_E, temp_NS, muRT)
    % Check condensed species
    aux = [];
    FLAG = false;
    if any(isnan(x))
        aux = true;
    else
        for i=length(temp_ind_swt_0):-1:1
            temp_NE = length(temp_ind_E);
            dG_dn = muRT(temp_ind_swt_0(i)) - dot(x(end-temp_NE:end-1), A0(temp_ind_swt_0(i), temp_ind_E));
            if dG_dn < 0
                aux = [aux, i];
            end
        end
    end
    if any(aux)
        FLAG = true;
    end
    if ~isempty(temp_ind_swt_0)
        temp_ind_swt = temp_ind_swt_0(aux);
    else
        temp_ind_swt = [];
    end
    temp_ind = [temp_ind_nswt, temp_ind_swt];
end

function ind_A = find_ind_Matrix(A, bool)
    ls = find(bool>0);
    ind_A = [];
    i = 1;
    for ind = ls
        ind_A = [ind_A, find(A(:, ind) > 0)'];
        i = i + 1;
    end
end

function ind_A0_E0 = remove_elements(NatomE, A0, tol)
    % Find zero sum elements
    bool_E0 = NatomE <= tol;
    ind_A0_E0 = find_ind_Matrix(A0, bool_E0);
end

function [temp_ind_nswt, temp_ind_swt, temp_ind_cryogenic, temp_ind_E, temp_NE] = temp_values(S, NatomE, tol)
    % List of indices with nonzero values and lengths
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_ind_cryogenic = S.ind_cryogenic;
    temp_NE = length(temp_ind_E);
end

function [ls1, ls2] = remove_item(N0, n, ind, ls1, ls2, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            if N0(ind(i), 2)
                ls1(ls1==ind(i)) = [];
            else
                ls2(ls2==ind(i)) = [];
            end
        end
    end
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [temp_ind, temp_ind_swt, temp_NS] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic)
    temp_ind = temp_ind(~ismember(temp_ind, temp_ind_cryogenic));
    temp_ind_swt = temp_ind_swt(~ismember(temp_ind_swt, temp_ind_cryogenic));
    temp_NS = length(temp_ind);
end

function lambda = relax_factor(NP, n, n_log_new, DeltaNP, temp_NG, SIZE)
    % Compute relaxation factor
    bool = log(n/NP) <= -SIZE & n_log_new >= 0;
    lambda = 2./max(5*abs(DeltaNP), abs(n_log_new));          
    lambda(bool(1:temp_NG)) = abs(-log(n(bool(1:temp_NG))/NP) - 9.2103404 ./ (n_log_new(bool(1:temp_NG)) - DeltaNP));
    lambda = min(1, min(lambda));
end

function [N0, NP] = apply_antilog(N0, NP_log, temp_ind_nswt)
    N0(temp_ind_nswt, 1) = exp(N0(temp_ind_nswt, 1));
    NP = exp(NP_log);
end

function STOP = compute_STOP(NP_0, NP, DeltaNP, N0, DeltaN0, temp_NG, A0, NatomE_tol, max_NatomE, tolE)
    DeltaN1 = N0 .* abs(DeltaN0) / NP;
    DeltaN1(temp_NG+1:end) = abs(DeltaN0(temp_NG+1:end)) / NP;
    DeltaN2 = NP_0 * abs(DeltaNP) / NP;
    Deltab = max(abs(NatomE_tol - sum(N0(:, 1) .* A0))) * max_NatomE;
    if Deltab < tolE
        Deltab = 0;
    end
    STOP = max(max(max(DeltaN1), DeltaN2), Deltab);
end

function A11 = update_matrix_A11(A0_T, N0, temp_ind_nswt, temp_ind_E, temp_NE)
    % Compute submatrix A11
    for k = temp_NE:-1:1
        A11(:, k) = sum(A0_T(temp_ind_E(k), temp_ind_nswt) .* A0_T(temp_ind_E, temp_ind_nswt) .* N0(temp_ind_nswt), 2);
    end
end

function A12 = update_matrix_A12(A0_T, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E)
    % Compute submatrix A12
    A12_1 = A0_T(temp_ind_E, temp_ind_swt);
    A12_2 = sum(A0_T(temp_ind_E, temp_ind_nswt) .* N0(temp_ind_nswt), 2);
    A12 = [A12_1, A12_2];
end

function A21 = update_matrix_A21(A0_T, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E)
    % Compute submaxtrix A21
    A21_1 = A0_T(temp_ind_E, temp_ind_swt)';
    A21_2 = sum(A0_T(temp_ind_E, temp_ind_nswt) .* N0(temp_ind_nswt), 2)';
    A21 = [A21_1; A21_2];
end

function A22 = update_matrix_A22(N0, NP, temp_ind_nswt, temp_NS, temp_NG)
    % Compute submatrix A22
    A22 = zeros(temp_NS - temp_NG + 1);
    A22(end, end) = sum(N0(temp_ind_nswt, 1) - NP);
end

function A = update_matrix_A(A0_T, N0, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS, temp_NE)
    % Compute matrix A
    A11 = update_matrix_A11(A0_T, N0, temp_ind_nswt, temp_ind_E, temp_NE);
    A12 = update_matrix_A12(A0_T, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E);
    A21 = update_matrix_A21(A0_T, N0, temp_ind_nswt, temp_ind_swt, temp_ind_E);
    A22 = update_matrix_A22(N0, NP, temp_ind_nswt, temp_NS, temp_NG);
    A = [A11, A12; A21, A22];
end

function b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, muRT) 
    % Compute vector b
    bi_0 = (NatomE(temp_ind_E) - N0(temp_ind, 1)' * A0(temp_ind, temp_ind_E) + sum(A0(temp_ind_nswt, temp_ind_E) .* N0(temp_ind_nswt, 1) .* muRT(temp_ind_nswt)))';    
    NP_0 = NP - sum(N0(temp_ind_nswt, 1)) + sum(N0(temp_ind_nswt, 1) .* muRT(temp_ind_nswt));
    b = [bi_0; muRT(temp_ind_swt); NP_0];
end

function Delta_ln_nj = update_Delta_ln_nj(A0, pi_i, Delta_NP, muRT, temp_ind_nswt, temp_ind_E)
    % Compute correction moles of gases
    Delta_ln_nj = sum(A0(temp_ind_nswt, temp_ind_E)' .* pi_i)' + Delta_NP - muRT(temp_ind_nswt);
end
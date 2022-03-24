function [N0, STOP] = equilibrium(self, pP, TP, strR)
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
    NatomE = strR.NatomE;
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
    temp_NS0 = temp_NS + 1;
    
    temp_ind_swt_0 = temp_ind_swt;
    temp_ind_swt = [];
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NS = length(temp_ind);
    % Initialize species vector N0
    N0(temp_ind_nswt, 1) = NP_0/temp_NG;
    % Standard Gibbs free energy
    g0 = set_g0(S.LS, TP, self.DB);
    % Dimensionless Chemical potential
    muRT = g0/R0TP;
    % Construction of part of matrix A
    [A1, temp_NS0] = update_matrix_A1(A0, [], temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
    A22 = zeros(temp_NE + 1);
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
            temp_NS0 = temp_NS + 1;
            % Initialize species vector N0 
            N0(:, 1) = NP_0/temp_NS;
            % Construction of part of matrix A
            [A1, temp_NS0] = update_matrix_A1(A0, [], temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
            A22 = zeros(temp_NE + 1);
            A0_T = A0';
        end
        STOP = 1;
        temp_NS = length(temp_ind);
        [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
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
            A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG);
            % Construction of vector b            
            b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, muRT);
            % Solve of the linear system A*x = b
            x = A\b;
            % Calculate correction factor
            lambda = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), x(end), temp_NG, SIZE);
            % Apply correction
            N0(temp_ind_nswt, 1) = log(N0(temp_ind_nswt, 1)) + lambda * x(1:temp_NG);
            N0(temp_ind_swt, 1) = N0(temp_ind_swt, 1) + lambda * x(temp_NG+1:temp_NS);
            NP_log = log(NP) + lambda * x(end);
            % Apply antilog
            [N0, NP] = apply_antilog(N0, NP_log, temp_ind_nswt);
            % Update temp values in order to remove species with moles < tolerance
            [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, NP, SIZE);
            % Update matrix A
            [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
            % Compute STOP criteria
            STOP = compute_STOP(NP_0, NP, x(end), N0(temp_ind, 1), x(1:temp_NS), temp_NG);
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

function [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A1
    A11 = eye(temp_NS);
    A11(temp_NG+1:end, temp_NG+1:end) = 0;
    A12 = -[A0(temp_ind, temp_ind_E), [ones(temp_NG, 1); zeros(temp_NS-temp_NG, 1)]];
    A1 = [A11, A12];
    temp_NS0 = temp_NS;
end

function [temp_ind, temp_ind_swt, temp_NS] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic)
    temp_ind = temp_ind(~ismember(temp_ind, temp_ind_cryogenic));
    temp_ind_swt = temp_ind(~ismember(temp_ind_swt, temp_ind_cryogenic));
    temp_NS = length(temp_ind);
end

function A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG)
    % Update stoichiometric submatrix A2
    A20 = N0(temp_ind, 1)';
    A20(temp_NG+1:end) = 0;
    A21 = [[N0(temp_ind_nswt, 1)' .* A0_T(temp_ind_E, temp_ind_nswt), A0_T(temp_ind_E, temp_ind_swt)]; A20];
    A22(end, end) = -NP;
    A2 = [A21, A22];
end

function A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG)
    % Update stoichiometric matrix A
    A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG);
    A = [A1; A2];
end

function b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, muRT)
    % Update coefficient vector b
    bi_0 = (NatomE(temp_ind_E) - N0(temp_ind, 1)' * A0(temp_ind, temp_ind_E))';
    NP_0 = NP - sum(N0(temp_ind_nswt, 1));
    b = [-muRT(temp_ind); bi_0; NP_0];
end

function relax = relax_factor(NP, n, n_log_new, DeltaNP, temp_NG, SIZE)
    % Compute relaxation factor
    bool = false(1, length(n));
    bool(1:temp_NG) = log(n(1:temp_NG))/log(NP) <= -SIZE & n_log_new(1:temp_NG) >= 0;
    lambda = ones(length(n), 1);
    lambda(~bool) = min(2./max(5*abs(DeltaNP), abs(n_log_new(~bool))), exp(2));          
    lambda(bool) = abs(-log(n(bool)/NP) - 9.2103404 ./ (n_log_new(bool) - DeltaNP));
    relax = min(1, min(lambda));  
end

function [N0, NP] = apply_antilog(N0, NP_log, temp_ind_nswt)
    N0(temp_ind_nswt, 1) = exp(N0(temp_ind_nswt, 1));
    NP = exp(NP_log);
end

function DeltaN = compute_STOP(NP_0, NP, DeltaNP, N0, DeltaN0, temp_NG)
    DeltaN1 = max(max(N0 .* abs(DeltaN0) / NP));
    DeltaN1(temp_NG+1:end) = abs(DeltaN0(temp_NG+1:end)) / NP;
    DeltaN3 = NP_0 * abs(DeltaNP) / NP;
    % Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
    DeltaN = max(DeltaN1, DeltaN3);
end
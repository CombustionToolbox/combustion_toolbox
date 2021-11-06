function [N0, STOP, STOP_ions] = equilibrium_ions(self, pP, TP, strR)
% Generalized Gibbs minimization method

% Abbreviations ---------------------
E = self.E;
S = self.S;
C = self.C;
TN = self.TN;
% -----------------------------------
% Set List of Species to List of Products
C.A0.value = C.A0.value(self.Misc.index_LS_original, :);
C.M0.value = C.M0.value(self.Misc.index_LS_original, :);
C.N0.value = C.N0.value(self.Misc.index_LS_original, :);
% -----------------------------------
N0 = C.N0.value;
A0 = C.A0.value;
R0TP = C.R0 * TP; % [J/(mol)]
% Initialization
% NatomE = N_CC(:,1)' * A0;
NatomE = strR.NatomE;
NP_0 = 0.1;
NP = NP_0;

it = 0;
itMax = 50 + round(S.NS/2);
SIZE = -log(TN.tolN);
STOP = 1.;
flag_ions_first = true;
% Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
% for the sum of elements whose value is <= tolN
temp_ind_ions = contains(S.LS, 'minus') | contains(S.LS, 'plus');
aux = NatomE;
if any(temp_ind_ions)
    NatomE(E.ind_E) = 1; % temporal fictitious value
end
ind_A0_E0 = remove_elements(NatomE, A0, TN.tolN);
NatomE = aux;
% List of indices with nonzero values
[temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_E, temp_NE] = temp_values(E.ind_E, S, NatomE, TN.tolN);
% Update temp values
[temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NS, ~, ~, ~] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, temp_ind_E, E, temp_ind_ions, [], TN.tolN, SIZE, flag_ions_first);
temp_NS0 = temp_NS + 1;
% Initialize species vector N0 
N0(temp_ind, 1) = 0.1/temp_NS;
% Dimensionless Standard Gibbs free energy 
g0 = set_g0(S.LS, TP, self.DB);
G0RT = g0/R0TP;
% Construction of part of matrix A (complete)
[A1, temp_NS0] = update_matrix_A1(A0, [], temp_NS, temp_NS0, temp_ind, temp_ind_E);
A22 = zeros(temp_NE + 1);
A0_T = A0';

while STOP > TN.tolN && it < itMax 
    it = it + 1;
    % Gibbs free energy
    G0RT(temp_ind_nswt) =  g0(temp_ind_nswt) / R0TP + log(N0(temp_ind_nswt, 1) / NP) + log(pP);
    % Construction of matrix A
    A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E);
    % Construction of vector b            
    b = update_vector_b(A0, N0, NP, NatomE, E.ind_E, temp_ind_ions, temp_ind, temp_ind_E, temp_ind_nswt, G0RT);
    % Solve of the linear system A*x = b
    x = A\b;
    % Compute correction factor
    lambda = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), x(end), SIZE);
    % Apply correction
    N0(temp_ind, 1) = log(N0(temp_ind, 1)) + lambda * x(1:temp_NS);
    NP_log = log(NP) + lambda * x(end);
    % Apply antilog
    [N0, NP] = apply_antilog(N0, NP_log, temp_ind);
    % Update temp values in order to remove species with moles < tolerance
    [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NS, temp_ind_E, A22, flag_ions_first] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, E, temp_ind_ions, A22, NP, SIZE, flag_ions_first);
    % Update matrix A
    [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NS, temp_NS0, temp_ind, temp_ind_E);
    % Compute STOP criteria
    STOP = compute_STOP(NP_0, NP, x(end), N0(temp_ind, 1), x(1:temp_NS));
end
% Check convergence of charge balance (ionized species)
[N0, STOP_ions] = check_convergence_ions(N0, A0, E.ind_E, temp_ind_nswt, temp_ind_ions, TN.tolN, TN.tol_pi_e);
% N0(N0(:, 1) < TN.tolN, 1) = 0;
end
% NESTED FUNCTIONS
function g0 = set_g0(ls, TP, DB)
    for i=length(ls):-1:1
        species = ls{i};
        g0(i, 1) = species_g0(species, TP, DB) * 1e3;
    end
end

function ind_A = find_ind_Matrix(A, bool)
    ls = find(bool>0);
%     ind_A = zeros(1, length(ls));
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

function [temp_ind_nswt, temp_ind_swt, temp_ind_ions, temp_ind_E, temp_NE] = temp_values(ind_E, S, NatomE, tol)
    % List of indices with nonzero values and lengths
    flag_ions = contains(S.LS, 'minus') | contains(S.LS, 'plus');
    if any(flag_ions)
        NatomE(ind_E) = 1; % temporal fictitious value
    end
    temp_ind_E    = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt  = S.ind_swt;
    temp_ind_ions = S.ind_nswt(flag_ions);
    temp_NE       = length(temp_ind_E);
end

function [temp_ind_swt, temp_ind_nswt, temp_ind_E, temp_ind_ions, A22, flag_ions_first] = remove_item(N0, n, ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, E, temp_ind_ions, A22, NP, SIZE, flag_ions_first)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            if N0(ind(i), 2)
                temp_ind_swt(temp_ind_swt==ind(i)) = [];
            else
                temp_ind_nswt(temp_ind_nswt==ind(i)) = [];
                temp_ind_ions(temp_ind_ions==ind(i)) = [];
                if isempty(temp_ind_ions) && flag_ions_first
                    % remove element E from matrix
                    temp_ind_E(E.ind_E) = [];
                    A22 = zeros(length(temp_ind_E) + 1);
                    flag_ions_first = false;
                end
            end
        end
    end
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_ions, temp_NS, temp_ind_E, A22, flag_ions_first] = update_temp(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, temp_ind_E, E, temp_ind_ions, A22, NP, SIZE, flag_ions_first)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt, temp_ind_E, temp_ind_ions, A22, flag_ions_first] = remove_item(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, temp_ind_E, E, temp_ind_ions, A22, NP, SIZE, flag_ions_first);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NS = length(temp_ind);
end
function [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NS, temp_NS0, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A1
    if temp_NS < temp_NS0
        A11 = eye(temp_NS);
        A12 = -[A0(temp_ind, temp_ind_E), ones(temp_NS, 1)];
        A1 = [A11, A12];
        temp_NS0 = temp_NS;
    end
end
function A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A2
    A21 = [N0(temp_ind, 1)' .* A0_T(temp_ind_E, temp_ind); N0(temp_ind, 1)'];
    A22(end, end) = -NP;
    A2 = [A21, A22];
end

function A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E)
    % Update stoichiometric matrix A
    A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_E);
    A = [A1; A2];
end

function b = update_vector_b(A0, N0, NP, NatomE, ind_E, temp_ind_ions, temp_ind, temp_ind_E, temp_ind_nswt, G0RT)
    % Update coefficient vector b
    bi_0 = (NatomE(temp_ind_E) - N0(temp_ind, 1)' * A0(temp_ind, temp_ind_E))';
    if any(temp_ind_ions)
        bi_0(temp_ind_E == ind_E) = 0;
    end
    NP_0 = NP - sum(N0(temp_ind_nswt, 1));
    b = [-G0RT(temp_ind); bi_0; NP_0];
end

function relax = relax_factor(NP, n, n_log_new, DeltaNP, SIZE)
    % Compute relaxation factor
    bool = log(n)/log(NP) <= -SIZE & n_log_new >= 0;
    lambda = ones(length(n), 1);
    lambda(bool) = abs(-log(n(bool)/NP) - 9.2103404 ./ (n_log_new(bool) - DeltaNP));
    lambda(~bool) = min(2./max(5*abs(DeltaNP), abs(n_log_new(~bool))), exp(2));          
    relax = min(1, min(lambda));  
end

function [relax, DeltaN3] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions)
    if any(temp_ind_ions)
        relax = -sum(A0(temp_ind_nswt, ind_E) .* N0(temp_ind_nswt, 1))/ ...
                 sum(A0(temp_ind_nswt, ind_E).^2 .* N0(temp_ind_nswt, 1));
        DeltaN3 = abs(sum(N0(temp_ind_nswt, 1) .* A0(temp_ind_nswt, ind_E)));
    else
        relax = [];
        DeltaN3 = 0;
    end 
end

function [N0, NP] = apply_antilog(N0, NP_log, temp_ind)
    N0(temp_ind, 1) = exp(N0(temp_ind, 1));
    NP = exp(NP_log);
end

function [DeltaN] = compute_STOP(NP_0, NP, DeltaNP, zip1, zip2)
    DeltaN1 = max(max(zip1 .* abs(zip2) / NP));
    DeltaN2 = NP_0 * abs(DeltaNP) / NP;
    DeltaN  = max(DeltaN1, DeltaN2);
end

function [N0, STOP] = check_convergence_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, TOL, TOL_pi)
    STOP = 0;
    if any(temp_ind_ions)
        [lambda_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
        if abs(lambda_ions) > TOL_pi
            [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, lambda_ions, TOL, TOL_pi);
        end
    end
end

function [N0, STOP] = recompute_ions(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions, lambda_ions, TOL, TOL_pi)
    % Parameters
    A0_ions = A0(temp_ind_ions, ind_E);
    STOP = 1;
    itmax = 30;
    it = 0;
    while STOP > TOL_pi && it < itmax
        it = it + 1;
        % Apply correction
        N0_ions =  log(N0(temp_ind_ions, 1)) + A0_ions * lambda_ions;
        % Apply antilog
        N0_ions = exp(N0_ions);
        % Assign values to original vector
        N0(temp_ind_ions, 1) = N0_ions;
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [lambda_ions, ~] = ions_factor(N0, A0, ind_E, temp_ind_nswt, temp_ind_ions);
        STOP = abs(lambda_ions);
    end   
    
    Xi_ions = N0(temp_ind_ions, 1) / sum(N0(:, 1));
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end 
end
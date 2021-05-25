function [N0, STOP] = Equilibrium(app, N_CC, phi, pP, TP, vR)
% Generalized Gibbs minimization method

% Abbreviations ---------------------
S = app.S;
C = app.C;
strThProp = app.strThProp;
% -----------------------------------

N0 = C.N0.value;
A0 = C.A0.value;
R0TP = C.R0 * TP; % [J/(mol)]
% Initialization
NatomE = N_CC(:,1)' * A0;
NP_0 = 0.1;
NP = NP_0;

it = 0;
% itMax = 500
itMax = 50 + round(S.NS/2);
SIZE = -log(C.tolN);
e = 0.;
STOP = 1.;

% Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
% for the sum of elements whose value is <= tolN
ind_A0_E0 = remove_elements(NatomE, A0, C.tolN);
% List of indices with nonzero values
[temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE] = temp_values(S, NatomE, C.tolN);
% Update temp values
[temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, C.tolN, SIZE);
% Initialize species vector N0 
N0(temp_ind, 1) = 0.1/temp_NS;
% Dimensionless Standard Gibbs free energy 
g0 = set_g0(S.LS, TP, strThProp);
G0RT = g0/R0TP;
% Construction of part of matrix A (complete)
A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E);
A22 = zeros(temp_NE + 1);
A0_T = A0';

while STOP > C.tolN && it < itMax
    it = it + 1;
    % Gibbs free energy
    G0RT(temp_ind_nswt) =  -(g0(temp_ind_nswt) / R0TP + log(N0(temp_ind_nswt, 1) / NP) + log(pP));
    % Construction of matrix A
    A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E);
    % Construction of vector b            
    b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, G0RT);
    % Solve of the linear system A*x = b
    x = A\b;
    % Calculate correction factor
    % update_SIZE(N0, A0, temp_ind, temp_ind_E, C.tolN)
    e = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), x(end), SIZE);
    % Apply correction
    N0_log = log(N0(temp_ind, 1)) + e * x(1:temp_NS);
    NP_log = log(NP) + e * x(end);
    % Apply antilog
    [N0, NP] = apply_antilog(N0, N0_log, NP_log, temp_ind);
    % Update temp values in order to remove species with moles < tolerance
    [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, NP, SIZE);
    % Update matrix A
    A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E);
    % Compute STOP criteria
    STOP = compute_STOP(NP_0, NP, x(end), N0(temp_ind, 1), x(1:temp_NS));
end
end
% NESTED FUNCTIONS
function g0 = set_g0(ls, TP, strThProp)
    for i=length(ls):-1:1
        species = ls{i};
        g0(i, 1) = species_g0(species, TP, strThProp)* 1e3;
    end
end
function ind_A = find_ind_Matrix(A, bool)
    ls = find(bool>0);
    ind_A = zeros(1, length(ls));
    i = 1;
    for ind = ls
        ind_A(i) = find(A(:, ind) > 0);
        i = i + 1;
    end
end
function ind_A0_E0 = remove_elements(NatomE, A0, tol)
    % Find zero sum elements
    bool_E0 = NatomE <= tol;
    ind_A0_E0 = find_ind_Matrix(A0, bool_E0);
end
function [temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE] = temp_values(S, NatomE, tol)
    % List of indices with nonzero values and lengths
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_NE = length(temp_ind_E);
end
function [ls1, ls2] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for n = zip1
        for ind = zip2
            if log(n/NP) < -SIZE
                if N0(ind, 2)
                    ls1(ind) = [];
                else
                    ls2(ind) = [];
                end
            end
        end
    end
end
function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS] = update_temp(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE);
    temp_ind = sort([temp_ind_nswt, temp_ind_swt]);
    temp_NS = length(temp_ind);
end
function A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A1
    A11 = eye(temp_NS);
    A12 = -[A0(temp_ind, temp_ind_E), ones(temp_NS, 1)];
    A1 = [A11, A12];
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

function b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, G0RT)
    % Update coefficient vector b
    for i=length(temp_ind_E):-1:1
        E = temp_ind_E(i);
        bi_0(i, 1) = NatomE(E) - N0(temp_ind, 1)' * A0(temp_ind, E);
    end
    NP_0 = NP - sum(N0(temp_ind_nswt, 1));
    b = [G0RT(temp_ind); bi_0; NP_0];
end

function relax = relax_factor(NP, zip1, zip2, DeltaNP, SIZE)
    % Compute relaxation factor
    e = [];
    for i=1:length(zip1)
        n = zip1(i);
        for j=1:length(zip2)
            n_log_new = zip2(j);
            if log(n)/log(NP) <= -SIZE && n_log_new >= 0.
                e = [e, abs(-log(n/NP) - 9.2103404 / (n_log_new - DeltaNP))];
            else
                e = [e, min(2/max(5*abs(DeltaNP), abs(n_log_new)), exp(2))];          
            end
        end
    end
    relax = min(1, min(e));  
end

function [N0, NP] = apply_antilog(N0, N0_log, NP_log, temp_ind)
    N0(temp_ind, :) = [exp(N0_log), N0(temp_ind, 2)];
    NP = exp(NP_log);
end

function DeltaN = compute_STOP(NP_0, NP, DeltaNP, zip1, zip2)
    for i=length(zip1):-1:1
        n = zip1(i);
        for j=length(zip2):-1:1
            n_log = zip2(i);
            DeltaN1(i,j) = n * abs(n_log) / NP;
        end
    end
    DeltaN1 = max(max(DeltaN1));
    DeltaN3 = NP_0 * abs(DeltaNP) / NP;
    % Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
    DeltaN = max(DeltaN1, DeltaN3);
end
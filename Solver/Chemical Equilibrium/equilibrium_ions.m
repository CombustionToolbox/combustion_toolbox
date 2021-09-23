function [N0, STOP] = equilibrium_ions(app, pP, TP, strR)
% Generalized Gibbs minimization method

% Abbreviations ---------------------
E = app.E;
S = app.S;
C = app.C;
TN = app.TN;
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
% itMax = 500;
itMax = 50 + round(S.NS/2);
SIZE = -log(TN.tolN);
STOP = 1.;
STOP_ions = 1.;
flag_ions_first = true;
% Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
% for the sum of elements whose value is <= tolN
flag_ions = contains(S.LS, 'minus') | contains(S.LS, 'plus');
aux = NatomE;
if any(flag_ions)
    NatomE(E.ind_E) = 1; % Fictitious value
end
ind_A0_E0 = remove_elements(NatomE, A0, TN.tolN);
NatomE = aux;
% List of indices with nonzero values
[temp_ind_nswt, temp_ind_swt, flag_ions, temp_ind_E, temp_NE] = temp_values(E.ind_E, S, NatomE, TN.tolN);
% Update temp values
[temp_ind, temp_ind_swt, temp_ind_nswt, flag_ions, temp_NS] = update_temp(N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, flag_ions, TN.tolN, SIZE);
temp_NS0 = temp_NS + 1;
% Initialize species vector N0 
N0(temp_ind, 1) = 0.1/temp_NS;
% Dimensionless Standard Gibbs free energy 
g0 = set_g0(S.LS, TP, app.strThProp);
G0RT = g0/R0TP;
% Construction of part of matrix A (complete)
[A1, temp_NS0] = update_matrix_A1(A0, [], temp_NS, temp_NS0, temp_ind, temp_ind_E);
A22 = zeros(temp_NE + 1);
A0_T = A0';

while (STOP > TN.tolN || STOP_ions > TN.tol_pi_e) && it < itMax 
    it = it + 1;
    % Gibbs free energy
    G0RT(temp_ind_nswt) =  g0(temp_ind_nswt) / R0TP + log(N0(temp_ind_nswt, 1) / NP) + log(pP);
    % Construction of matrix A
    A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_E);
    % Construction of vector b            
    b = update_vector_b(A0, N0, NP, NatomE, E.ind_E, flag_ions, temp_ind, temp_ind_E, temp_ind_nswt, G0RT);
    % Solve of the linear system A*x = b
    x = A\b;
    % Compute correction factor
    lambda = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), x(end), SIZE);
    % Compute and apply correction of the Lagrangian multiplier for ions divided by RT
    [lambda_ions, DeltaN3] = ions_factor(N0, A0, temp_ind_nswt, E.ind_E, flag_ions);
    % Apply correction
    if any(flag_ions) && flag_ions_first
        N0_wions =  log(N0(temp_ind_nswt(flag_ions), 1)) + A0(temp_ind_nswt(flag_ions), E.ind_E) * lambda_ions;
    end
    N0(temp_ind, 1) = log(N0(temp_ind, 1)) + lambda * x(1:temp_NS);
    if any(flag_ions)
        if abs(lambda_ions) > TN.tol_pi_e && flag_ions_first
            N0(temp_ind_nswt(flag_ions), 1) = N0_wions;
            flag_ions_first = false;
        end
    end
    NP_log = log(NP) + lambda * x(end);
    % Apply antilog
    [N0, NP] = apply_antilog(N0, NP_log, temp_ind);
    % Update temp values in order to remove species with moles < tolerance
    [temp_ind, temp_ind_swt, temp_ind_nswt, flag_ions, temp_NS] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, flag_ions, NP, SIZE);
    % Update matrix A
    [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NS, temp_NS0, temp_ind, temp_ind_E);
    % Compute STOP criteria
    [STOP, STOP_ions] = compute_STOP(NP_0, NP, x(end), N0(temp_ind, 1), x(1:temp_NS), lambda_ions, flag_ions, flag_ions_first, DeltaN3);
end
% N0(N0(:, 1) < TN.tolN, 1) = 0;
end
% NESTED FUNCTIONS
function g0 = set_g0(ls, TP, strThProp)
    for i=length(ls):-1:1
        species = ls{i};
        g0(i, 1) = species_g0(species, TP, strThProp) * 1e3;
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

function [temp_ind_nswt, temp_ind_swt, flag_ions, temp_ind_E, temp_NE] = temp_values(ind_E, S, NatomE, tol)
    % List of indices with nonzero values and lengths
    flag_ions = contains(S.LS, 'minus') | contains(S.LS, 'plus');
    if any(flag_ions)
        NatomE(ind_E) = 1; % Fictitious value
    end
    
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_NE = length(temp_ind_E);
end

function [temp_ind_swt, temp_ind_nswt, flag_ions] = remove_item(N0, n, ind, temp_ind_swt, temp_ind_nswt, flag_ions, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            if N0(ind(i), 2)
                temp_ind_swt(temp_ind_swt==ind(i)) = [];
            else
                temp_ind_nswt(temp_ind_nswt==ind(i)) = [];
                try
                    flag_ions(ind(i)) = [];
                catch
                    continue
                end
            end
        end
    end
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, flag_ions, temp_NS] = update_temp(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, flag_ions, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt, flag_ions] = remove_item(N0, zip1, zip2, temp_ind_swt, temp_ind_nswt, flag_ions, NP, SIZE);
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

function b = update_vector_b(A0, N0, NP, NatomE, ind_E, flag_ions, temp_ind, temp_ind_E, temp_ind_nswt, G0RT)
    % Update coefficient vector b
    bi_0 = (NatomE(temp_ind_E) - N0(temp_ind, 1)' * A0(temp_ind, temp_ind_E))';
    if any(flag_ions)
        bi_0(ind_E) = 0;
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

function [relax, DeltaN3] = ions_factor(N0, A0, temp_ind_nswt, ind_E, flag_ions)
    if any(flag_ions)
        relax = -sum(A0(temp_ind_nswt, ind_E)    .* N0(temp_ind_nswt, 1))/ ...
                 sum(A0(temp_ind_nswt, ind_E).^2 .* N0(temp_ind_nswt, 1));
        DeltaN3 = abs(sum(N0(temp_ind_nswt, 1) .* A0(temp_ind_nswt, ind_E)));
    else
        relax = [];
    end 
end

function [N0, NP] = apply_antilog(N0, NP_log, temp_ind)
    N0(temp_ind, 1) = exp(N0(temp_ind, 1));
    NP = exp(NP_log);
end

function [DeltaN, Delta_ions] = compute_STOP(NP_0, NP, DeltaNP, zip1, zip2, lambda_ions, flag_ions, flag_ions_first, DeltaN3)
    DeltaN1 = max(max(zip1 .* abs(zip2) / NP));
    DeltaN2 = NP_0 * abs(DeltaNP) / NP;
    DeltaN  = max(DeltaN1, max(DeltaN2, DeltaN3));
    if flag_ions_first
        Delta_ions = 0;
    elseif any(flag_ions)
        Delta_ions = abs(lambda_ions);
    else
        Delta_ions = 0;
    end
end
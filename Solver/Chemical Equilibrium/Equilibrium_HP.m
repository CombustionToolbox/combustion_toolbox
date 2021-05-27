function [N0, STOP] = Equilibrium_HP(app, pP, TP, strR)
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
% NatomE = N_CC(:,1)' * A0;
NatomE = strR.NatomE;
NP_0 = 0.1;
NP = NP_0;

it = 0;
itMax = 500;
% itMax = 50 + round(S.NS/2);
SIZE = -log(C.tolN);
e = 0.;
STOP = 1.;
STOP_HP = 1.;
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
cp0 = set_cp0(S.LS, TP, strThProp);
H0 = set_H0(S.LS, TP, strThProp);
g0 = set_g0(S.LS, TP, strThProp);
CP0R = cp0/C.R0;
H0RT = H0/R0TP;
h0 = strR.h*1e3;
% Construction of part of matrix A (complete)
A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E, H0RT);
A22 = zeros(temp_NE + 2);
A0_T = A0';

% while STOP > C.tolN && it < itMax && STOP_HP > 1e-4
while STOP_HP > 1e-4
    it = it + 1;
    % Enthalpy
    h0RT = h0 / (C.R0 * TP);
    H0RT = H0 / (C.R0 * TP);
    % Gibbs free energy
    G0RT = g0 / (C.R0 * TP);
    G0RT(temp_ind_nswt) =  g0(temp_ind_nswt) / (C.R0 * TP) + log(N0(temp_ind_nswt, 1) / NP) + log(pP);
    % Construction of matrix A
    A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, CP0R, H0RT);
    % Construction of vector b            
    b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, h0RT, H0RT, G0RT);
    % Solve of the linear system A*x = b
    x = A\b;
    % Calculate correction factor
    % update_SIZE(N0, A0, temp_ind, temp_ind_E, C.tolN)
    e = relax_factor(NP, N0(temp_ind, 1), x(1:temp_NS), x(end-1), x(end), SIZE);
    % Apply correction
    N0_log = [log(N0(temp_ind_nswt, 1)) + e * x(1:temp_NG), N0(temp_ind_swt, 1) + e * x(temp_NG:temp_NS)];
    NP_log = log(NP) + e * x(end-1);
    T_log = log(TP) + e * x(end);
    % Apply antilog
    [N0, NP, TP] = apply_antilog(N0, N0_log, NP_log, T_log, temp_ind);
    % Update temp values in order to remove species with moles < tolerance
    [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(N0, N0(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, NP, SIZE);
    % Update matrix A1
    A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E, H0RT);
    % Compute STOP criteria
    STOP = compute_STOP(NP_0, NP, x(end-1), N0(temp_ind, 1), x(1:temp_NS));
    STOP_HP = abs(T_log);
end
end
% NESTED FUNCTIONS
function cp0 = set_cp0(ls, TP, strThProp)
    for i=length(ls):-1:1
        cp0(i, 1) = species_cP(ls{i}, TP, strThProp);
    end
end

function H0 = set_H0(ls, TP, strThProp)
    for i=length(ls):-1:1
%         h0(i, 1) = species_DhT(ls{i}, TP, strThProp)* 1e3;
        H0(i, 1) = strThProp.(ls{i}).hf * 1e3;
    end
end

function g0 = set_g0(ls, TP, strThProp)
    for i=length(ls):-1:1
        g0(i, 1) = species_g0(ls{i}, TP, strThProp)* 1e3;
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

function [temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE] = temp_values(S, NatomE, tol)
    % List of indices with nonzero values and lengths
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
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
%     bool = log(zip1./NP) < -SIZE;
%     bool1 = bool & N0(zip2, 2);
%     bool2 = bool & ~N0(zip2, 2);
%     ind1 = logical(sum(ls1(:)==zip2(bool1), 2)); if ~ind1, ind1=[]; end
%     ind2 = logical(sum(ls2(:)==zip2(bool2), 2)); if ~ind2, ind2=[]; end
%     ls1(ind1) = [];
%     ls2(ind2) = [];
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end
function A1 = update_matrix_A1(A0, temp_NS, temp_ind, temp_ind_E, H0RT)
    % Update stoichiometric submatrix A1
    A11 = eye(temp_NS);
    A12 = -[A0(temp_ind, temp_ind_E), ones(temp_NS, 1)];
    A13 = -H0RT(temp_ind);
    A1 = [A11, A12, A13];
end
function A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, CP0R, H0RT)
    % Update stoichiometric submatrix A2
    A20 = [N0(temp_ind_nswt, 1) .* H0RT(temp_ind_nswt); H0RT(temp_ind_swt)]';
    A21 = [N0(temp_ind, 1)' .* A0_T(temp_ind_E, temp_ind); N0(temp_ind, 1)'; A20];
    A22(end-1, end-1) = -NP;
    A22(end, end) = sum(N0(temp_ind, 1) .* CP0R(temp_ind));
    A2 = [A21, A22];
end

function A = update_matrix_A(A0_T, A1, A22, N0, NP, temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, CP0R, H0RT)
    % Update stoichiometric matrix A
    A2 = update_matrix_A2(A0_T, A22, N0, NP, temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_E, CP0R, H0RT);
    A = [A1; A2];
end

function b = update_vector_b(A0, N0, NP, NatomE, temp_ind, temp_ind_E, temp_ind_nswt, h0RT, H0RT, G0RT)
    % Update coefficient vector b
    bi_0 = (NatomE(temp_ind_E) - N0(temp_ind, 1)' * A0(temp_ind, temp_ind_E))';
    NP_0 = NP - sum(N0(temp_ind_nswt, 1));
    b = [-G0RT(temp_ind); bi_0; NP_0; h0RT-sum(N0(temp_ind, 1) .* H0RT(temp_ind))];
end

function relax = relax_factor(NP, n, n_log_new, DeltaNP, DeltaT, SIZE)
    % Compute relaxation factor
    bool = log(n)/log(NP) <= -SIZE & n_log_new >= 0;
    e = ones(length(n), 1);
    e(bool) = abs(-log(n(bool)/NP) - 9.2103404 ./ (n_log_new(bool) - DeltaNP));
    e(~bool) = min(2./max(max(5*abs(DeltaT), 5*abs(DeltaNP)), abs(n_log_new(~bool))), exp(2));          
    relax = min(1, min(e));  
end

function [N0, NP, TP] = apply_antilog(N0, N0_log, NP_log, TP_log, temp_ind)
    N0(temp_ind_nswt, :) = [exp(N0_log(temp_ind_nswt)), N0(temp_ind_nswt, 2)];
    NP = exp(NP_log);
    TP = exp(TP_log);
end

function DeltaN = compute_STOP(NP_0, NP, DeltaNP, zip1, zip2)
    DeltaN1 = max(max(zip1 .* abs(zip2) / NP));
    DeltaN3 = NP_0 * abs(DeltaNP) / NP;
    % Deltab = [abs(bi - sum(N0[:, 0] * A0[:, i])) for i, bi in enumerate(x[S.NS:-1]) if bi > 1e-6]
    DeltaN = max(DeltaN1, DeltaN3);
end
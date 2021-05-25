function [N_IC, DeltaNP] = Equilibrium(app, N_CC, phi, pP, TP, vR)
% Generalized Gibbs minimization method

% Abbreviations ---------------------
S = app.S;
C = app.C;
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
[temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE, temp_NS, temp_ind_remove] = temp_values(S, NatomE, C.tolN);
% Update temp values
[temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS] = update_temp(temp_ind, temp_NS, N0, N0(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, C.tolN, SIZE);
end
% NESTED FUNCTIONS
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
function [temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NE, temp_NS, temp_ind_remove] = temp_values(S, NatomE, tol)
    % List of indices with nonzero values and lengths
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_ind = temp_ind_nswt + temp_ind_swt;
    temp_NE = length(temp_ind_E);
    temp_NS = S.NS;
    temp_ind_remove = [];
end
function [N0, ls0, ls1, ls2] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    ls0 = [];
    for n = zip1
        for ind = zip2
            if log(n/NP) < -SIZE
                ls0 = [ls0, ind];
                N0(ind, 1) = 0.;
                if N0(ind, 2)
                    ls1(ind) = [];
                else
                    ls2(ind) = [];
                end
            end
        end
    end
end
function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NS] = update_temp(temp_ind, temp_NS, N0, zip1, zip2, ls1, ls2, NP, SIZE)
    % Update temp items
    [N0, temp_ind_remove, temp_ind_swt, temp_ind_nswt] = remove_item(N0, zip1, zip2, ls1, ls2, NP, SIZE);
    if temp_ind_remove
        temp_ind = list(set(temp_ind) - set(temp_ind_remove));
        temp_NS = len(temp_ind);
    end
end
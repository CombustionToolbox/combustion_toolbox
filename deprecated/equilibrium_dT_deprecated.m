function [dNi_T, dN_T] = equilibrium_dT_deprecated(self, moles, T, mix1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_nswt, temp_ind_swt, temp_ind_E)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     self (struct):   Data of the mixture, conditions, and databases
    %     moles (float):   Equilibrium composition [moles]
    %     T (float):       Temperature [K]
    %     mix1 (struct):   Properties of the initial mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     - dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     - dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature

    % Abbreviations ---------------------
    E = self.E;
    S = self.S;
    C = self.C;
    TN = self.TN;
    % -----------------------------------
    A0 = C.A0.value;
    R0TP = C.R0 * T; % [J/mol]
    % Initialization
    NatomE = mix1.NatomE;
    NP = sum(moles(:, 1));
    dNi_T = zeros(length(moles), 1);
    SIZE = -log(TN.tolN);
    
    temp_NG = length(temp_ind_nswt);
    temp_NE = length(temp_ind_E);
%     % Find indeces of the species/elements that we have to remove from the stoichiometric matrix A0
%     % for the sum of elements whose value is <= tolN
%     temp_ind_ions = contains(S.LS, 'minus') | contains(S.LS, 'plus'); %S.ind_ions;
%     if any(temp_ind_ions)
%         pass_ions = moles(temp_ind_ions) > TN.tolN;
%         temp_ind_ions = temp_ind_ions(pass_ions);
%         aux = NatomE;
%         NatomE(E.ind_E) = 1; % temporal fictitious value
%     end
%     ind_A0_E0 = remove_elements(NatomE, A0, TN.tolN);
%     NatomE = aux;
%     % List of indices with nonzero values
%     [temp_ind_nswt, temp_ind_swt, temp_ind_cryogenic, temp_ind_E, temp_NE] = temp_values(S, NatomE, TN.tolN);
%     % Update temp values
%     [temp_ind, temp_ind_swt, temp_ind_nswt, ~, ~] = update_temp(moles, moles(ind_A0_E0, 1), ind_A0_E0, temp_ind_swt, temp_ind_nswt, temp_ind_cryogenic, NP, SIZE);
%     [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(moles, moles(temp_ind, 1), temp_ind, temp_ind_swt, temp_ind_nswt, temp_ind_cryogenic, NP, SIZE);
    temp_NS0 = temp_NS + 1;
    % Dimensionless Standard-state enthalpy [J/mol]
    h0 = set_h0(S.LS, T, self.DB);
    H0RT =  h0 / R0TP;
    % Construction of part of matrix A (complete)
    [A1, ~] = update_matrix_A1(A0, [], temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E);
    A22 = zeros(temp_NE + 1);
    A0_T = A0';
    % Construction of matrix A
    A = update_matrix_A(A0_T, A1, A22, moles, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS);
    % Construction of vector b            
    b = update_vector_b(temp_ind, temp_ind_E, H0RT);
    % Solve of the linear system A*x = b
    x = A\b;
    dNi_T(temp_ind) = x(1:temp_NS);
    dN_T = x(end);
end

% SUB-PASS FUNCTIONS
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

function [temp_ind_nswt, temp_ind_swt, temp_ind_cryogenic, temp_ind_E, temp_NE] = temp_values(S, NatomE, tol)
    % List of indices with nonzero values and lengths
    temp_ind_E = find(NatomE > tol);
    temp_ind_nswt = S.ind_nswt;
    temp_ind_swt = S.ind_swt;
    temp_ind_cryogenic = S.ind_cryogenic;
    temp_NE = length(temp_ind_E);
end

function [ls1, ls2] = remove_item(moles, n, ind, ls1, ls2, NP, SIZE)
    % Remove species from the computed indeces list of gaseous and condensed
    % species and append the indeces of species that we have to remove
    for i=1:length(n)
        if log(n(i)/NP) < -SIZE
            if moles(ind(i), 2)
                ls1(ls1==ind(i)) = [];
            else
                ls2(ls2==ind(i)) = [];
            end
        end
    end
end

function [temp_ind, temp_ind_swt, temp_ind_nswt, temp_NG, temp_NS] = update_temp(moles, zip1, zip2, ls1, ls2, temp_ind_cryogenic, NP, SIZE)
    % Update temp items
    [temp_ind_swt, temp_ind_nswt] = remove_item(moles, zip1, zip2, ls1, ls2, NP, SIZE);
    temp_ind = [temp_ind_nswt, temp_ind_swt];
    [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic);
    temp_NG = length(temp_ind_nswt);
    temp_NS = length(temp_ind);
end

function [A1, temp_NS0] = update_matrix_A1(A0, A1, temp_NG, temp_NS, temp_NS0, temp_ind, temp_ind_E)
    % Update stoichiometric submatrix A1
    if temp_NS < temp_NS0
        A11 = eye(temp_NS);
        A11(temp_NG+1:end, temp_NG+1:end) = 0;
        A12 = -[A0(temp_ind, temp_ind_E), [ones(temp_NG, 1); zeros(temp_NS-temp_NG, 1)]];
        A1 = [A11, A12];
        temp_NS0 = temp_NS;
    end
end

function [temp_ind, temp_ind_swt] = check_cryogenic(temp_ind, temp_ind_swt, temp_ind_cryogenic)
    temp_ind = temp_ind(~ismember(temp_ind, temp_ind_cryogenic));
    temp_ind_swt = temp_ind_swt(~ismember(temp_ind_swt, temp_ind_cryogenic));
end

function A2 = update_matrix_A2(A0_T, A22, moles, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS)
    % Update stoichiometric submatrix A2
    A21_1 = [moles(temp_ind_nswt, 1)' .* A0_T(temp_ind_E, temp_ind_nswt); moles(temp_ind_nswt, 1)'];
    A21_2 = [A0_T(temp_ind_E, temp_ind_swt); zeros(1, temp_NS-temp_NG)];
    A21 = [A21_1, A21_2];
    A22(end, end) = -NP;
    A2 = [A21, A22];
end

function A = update_matrix_A(A0_T, A1, A22, moles, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS)
    % Update stoichiometric matrix A
    A2 = update_matrix_A2(A0_T, A22, moles, NP, temp_ind_nswt, temp_ind_swt, temp_ind_E, temp_NG, temp_NS);
    A = [A1; A2];
end

function b = update_vector_b(temp_ind, temp_ind_E, H0RT)
    % Update coefficient vector b
    b = [H0RT(temp_ind); zeros(length(temp_ind_E) + 1, 1)];
end
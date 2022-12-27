function [dNi_T, dN_T, A] = equilibrium_dT_reduced(self, N0, T, A0, NG, NS, NE, ind_nswt, ind_swt, ind_elem)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     N0 (float): Equilibrium composition [moles]
    %     T (float): Temperature [K]
    %     A0 (float): Stoichiometric matrix
    %     NG (float): Temporal total number of gaseous species
    %     NS (float): Temporal total number of species
    %     NE (float): Temporal total number of elements
    %     ind (float): Temporal index of species in the final mixture
    %     ind_nswt (float): Temporal index of gaseous species in the final mixture
    %     ind_swt (float): Temporal index of condensed species in the final mixture
    %     ind_elem (float): Temporal index of elements in the final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature
    %     * A (float): Matrix A to solve the linear system A*x = b

    % Definitions
    R0TP = self.C.R0 * T; % [J/mol]
    % Initialization
    dNi_T = zeros(length(N0), 1);
    % Dimensionless Standard-state enthalpy [J/mol]
    h0 = set_h0(self.S.LS, T, self.DB);
    H0RT = h0 / R0TP;
    % Construction of part of matrix A
    A22 = zeros(NS - NG + 1);
    A0 = A0(:, ind_elem);
    A0_T = A0';
    % Construction of matrix A
    A = update_matrix_A(A0_T, A22, N0, ind_nswt, ind_swt, NE);
    % Construction of vector b
    b = update_vector_b(A0, N0, ind_nswt, ind_swt, H0RT);
    % Solve of the linear system A*x = b
    x = A \ b;
    % Extract solution
    dpii_T = x(1:NE);
    dNi_T(ind_swt) = x(NE+1:end-1);
    dN_T = x(end);
    % Compute remainder dNi_T (gas)
    dNi_T(ind_nswt) = H0RT(ind_nswt) + A0(ind_nswt, :) * dpii_T + dN_T;
end

% SUB-PASS FUNCTIONS
function A11 = update_matrix_A11(A0_T, N0, ind_nswt, NE)
    % Compute submatrix A11
    for k = NE:-1:1
        A11(:, k) = sum(A0_T(k, ind_nswt) .* A0_T(:, ind_nswt) .* N0(ind_nswt), 2);
    end
end

function A12 = update_matrix_A12(A0_T, N0, ind_nswt, ind_swt)
    % Compute submatrix A12
    A12_1 = A0_T(:, ind_swt);
    A12_2 = sum(A0_T(:, ind_nswt) .* N0(ind_nswt), 2);
    A12 = [A12_1, A12_2];
end

function A = update_matrix_A(A0_T, A22, N0, ind_nswt, ind_swt, NE)
    % Compute matrix A
    A11 = update_matrix_A11(A0_T, N0, ind_nswt, NE);
    A12 = update_matrix_A12(A0_T, N0, ind_nswt, ind_swt);
    A = [A11, A12; A12', A22];
end

function b = update_vector_b(A0, N0, ind_nswt, ind_swt, H0RT)
    % Compute vector b
    b = -[(sum(A0(ind_nswt, :) .* N0(ind_nswt, 1) .* H0RT(ind_nswt)))';...
        H0RT(ind_swt); sum(N0(ind_nswt, 1) .* H0RT(ind_nswt))];
end
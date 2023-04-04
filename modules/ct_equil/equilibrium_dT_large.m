function [dNi_T, dN_T] = equilibrium_dT_large(self, moles, T, A0, NG, NS, NE, ind, ind_nswt, ind_swt, ind_E)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     moles (float): Equilibrium composition [moles]
    %     T (float): Temperature [K]
    %     A0 (float): Stoichiometric matrix
    %     NG (float): Temporal total number of gaseous species
    %     NS (float): Temporal total number of species
    %     NE (float): Temporal total number of elements
    %     ind (float): Temporal index of species in the final mixture
    %     ind_nswt (float): Temporal index of gaseous species in the final mixture
    %     ind_swt (float): Temporal index of condensed species in the final mixture
    %     ind_E (float): Temporal index of elements in the final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature
    
    % Definitions
    R0TP = self.C.R0 * T; % [J/mol]

    % Initialization
    NP = sum(moles(ind_nswt, 1));
    dNi_T = zeros(length(moles), 1);

    % Molar enthalpy [J/mol]
    h0 = set_h0(self.S.LS, T, self.DB);

    % Dimensionless enthalpy
    H0RT = h0 / R0TP;

    % Construction of part of matrix J
    J1 = update_matrix_J1(A0, NG, NS, ind, ind_E);
    J22 = zeros(NE + 1);
    A0_T = A0';

    % Construction of matrix J
    J = update_matrix_J(A0_T, J1, J22, moles, NP, ind_nswt, ind_swt, ind_E, NG, NS);
    
    % Construction of vector b
    b = update_vector_b(ind, ind_E, H0RT);
    
    % Solve of the linear system J*x = b
    x = J \ b;

    % Extract solution
    dNi_T(ind) = x(1:NS);
    dN_T = x(end);
end

% SUB-PASS FUNCTIONS
function J1 = update_matrix_J1(A0, NG, NS, ind, ind_E)
    % Update stoichiometric submatrix J1
    J11 = eye(NS);
    J11(NG + 1:end, NG + 1:end) = 0;
    J12 =- [A0(ind, ind_E), [ones(NG, 1); zeros(NS - NG, 1)]];
    J1 = [J11, J12];
end

function J2 = update_matrix_J2(A0_T, J22, moles, NP, ind_nswt, ind_swt, ind_E, NG, NS)
    % Update stoichiometric submatrix J2
    J21_1 = [moles(ind_nswt, 1)' .* A0_T(ind_E, ind_nswt); moles(ind_nswt, 1)'];
    J21_2 = [A0_T(ind_E, ind_swt); zeros(1, NS - NG)];
    J21 = [J21_1, J21_2];
    J22(end, end) = -NP;
    J2 = [J21, J22];
end

function J = update_matrix_J(A0_T, J1, J22, moles, NP, ind_nswt, ind_swt, ind_E, NG, NS)
    % Update stoichiometric matrix J
    J2 = update_matrix_J2(A0_T, J22, moles, NP, ind_nswt, ind_swt, ind_E, NG, NS);
    J = [J1; J2];
end

function b = update_vector_b(ind, ind_E, H0RT)
    % Update coefficient vector b
    b = [H0RT(ind); zeros(length(ind_E) + 1, 1)];
end

function [dNi_T, dN_T] = equilibrium_dT(J, N0, A0, NE, ind_nswt, ind_swt, ind_elem, H0RT)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N0 (float): Equilibrium composition [moles]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     ind (float): Temporal index of species in the final mixture
    %     ind_nswt (float): Temporal index of gaseous species in the final mixture
    %     ind_swt (float): Temporal index of condensed species in the final mixture
    %     ind_elem (float): Temporal index of elements in the final mixture
    %     H0RT (float): Dimensionless standard-state enthalpy
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature

    % Definitions
    A0 = A0(:, ind_elem);
    % Initialization
    dNi_T = zeros(length(N0), 1);
    % Construction of vector b
    b = update_vector_b(A0, N0, ind_nswt, ind_swt, H0RT);
    % Solve of the linear system J*x = b
    x = J \ b;
    % Extract solution
    dpii_T = x(1:NE);
    dNi_T(ind_swt) = x(NE+1:end-1);
    dN_T = x(end);
    % Compute remainder dNi_T (gas)
    dNi_T(ind_nswt) = H0RT(ind_nswt) + A0(ind_nswt, :) * dpii_T + dN_T;
end

% SUB-PASS FUNCTIONS
function b = update_vector_b(A0, N0, ind_nswt, ind_swt, H0RT)
    % Compute vector b
    b = -[(sum(A0(ind_nswt, :) .* N0(ind_nswt, 1) .* H0RT(ind_nswt)))';...
        H0RT(ind_swt); sum(N0(ind_nswt, 1) .* H0RT(ind_nswt))];
end
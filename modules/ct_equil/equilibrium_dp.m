function [dNi_p, dN_p] = equilibrium_dp(J, N0, A0, NE, ind_nswt, ind_swt, ind_elem)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to pressure from a given composition [moles] at equilibrium
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
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float):  Thermodynamic derivative of the moles of the mixture respect to pressure

    % Definitions
    A0 = A0(:, ind_elem);
    % Initialization
    dNi_p = zeros(length(N0), 1);
    % Construction of vector b
    b = update_vector_b(J, N0, ind_nswt);
    % Solve of the linear system J*x = b
    x = J \ b;
    % Extract solution
    dpii_p = x(1:NE);
    dNi_p(ind_swt) = x(NE+1:end-1);
    dN_p = x(end);
    % Compute remainder dNi_T (gas)
    dNi_p(ind_nswt) = -1 + A0(ind_nswt, :) * dpii_p + dN_p;
end

% SUB-PASS FUNCTIONS
function b = update_vector_b(J, N0, ind_nswt)
    % Compute vector b
    b = J(:, end);
    b(end) = sum(N0(ind_nswt, 1));
end
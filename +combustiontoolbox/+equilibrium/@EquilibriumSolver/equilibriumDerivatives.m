function [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N0, A0, NE, indexGas, indexCondensed, indexElements, H0RT)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature and pressure from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N0 (float): Equilibrium composition [moles]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
    %     indexElements (float): Temporal index of elements in the final mixture
    %     H0RT (float): Dimensionless enthalpy
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float):  Thermodynamic derivative of the moles of the mixture respect to pressure
    %
    % Example:
    %     [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N0, A0, NE, ind, indexGas, indexCondensed, indexElements, H0RT)
    
    % Equilibrium derivative respect temperature
    [dNi_T, dN_T] = equilibrium_dT(J, N0, A0, NE, indexGas, indexCondensed, indexElements, H0RT);

    % Equilibrium derivative respect pressure
    [dNi_p, dN_p] = equilibrium_dp(J, N0, A0, NE, indexGas, indexCondensed, indexElements);
end

function [dNi_T, dN_T] = equilibrium_dT(J, N0, A0, NE, indexGas, indexCondensed, indexElements, H0RT)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N0 (float): Equilibrium composition [moles]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
    %     indexElements (float): Temporal index of elements in the final mixture
    %     H0RT (float): Dimensionless enthalpy
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature
    %
    % Example:
    %     [dNi_T, dN_T] = equilibrium_dT(J, N0, A0, NE, ind, indexGas, indexCondensed, indexElements, H0RT)

    % Definitions
    A0 = A0(:, indexElements);

    % Initialization
    dNi_T = zeros(length(N0), 1);

    % Construction of vector b
    b = update_vector_b(A0, N0, indexGas, indexCondensed, H0RT);

    % Solve of the linear system J*x = b
    x = J \ b;

    % Extract solution
    dpii_T = x(1:NE);
    dNi_T(indexCondensed) = x(NE+1:end-1);
    dN_T = x(end);

    % Compute remainder dNi_T (gas)
    dNi_T(indexGas) = H0RT(indexGas) + A0(indexGas, :) * dpii_T + dN_T;

    % NESTED FUNCTION
    function b = update_vector_b(A0, N0, indexGas, indexCondensed, H0RT)
        % Compute vector b
        b = -[(sum(A0(indexGas, :) .* N0(indexGas, 1) .* H0RT(indexGas)))';...
              H0RT(indexCondensed); sum(N0(indexGas, 1) .* H0RT(indexGas))];
    end

end

function [dNi_p, dN_p] = equilibrium_dp(J, N0, A0, NE, indexGas, indexCondensed, indexElements)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to pressure from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N0 (float): Equilibrium composition [moles]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
    %     indexElements (float): Temporal index of elements in the final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float):  Thermodynamic derivative of the moles of the mixture respect to pressure
    %
    % Example:
    %     [dNi_p, dN_p] = equilibrium_dp(J, N0, A0, NE, ind, ind_nswt, ind_swt, ind_elem)

    % Definitions
    A0 = A0(:, indexElements);

    % Initialization
    dNi_p = zeros(length(N0), 1);

    % Construction of vector b
    b = update_vector_b(J, N0, indexGas);

    % Solve of the linear system J*x = b
    x = J \ b;

    % Extract solution
    dpii_p = x(1:NE);
    dNi_p(indexCondensed) = x(NE+1:end-1);
    dN_p = x(end);
    
    % Compute remainder dNi_T (gas)
    dNi_p(indexGas) = -1 + A0(indexGas, :) * dpii_p + dN_p;

    % NESTED FUNCTION
    function b = update_vector_b(J, N0, ind_nswt)
        % Compute vector b
        b = J(:, end);
        b(end) = sum(N0(ind_nswt, 1));
    end

end
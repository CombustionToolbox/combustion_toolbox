function [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N, A0, NE, indexGas, indexCondensed, H0RT)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature and pressure from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
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
    %     [dNi_T, dN_T, dNi_p, dN_p] = equilibriumDerivatives(J, N, A0, NE, ind, indexGas, indexCondensed, H0RT)
    
    % Equilibrium derivative respect temperature
    [dNi_T, dN_T] = equilibrium_dT(J, N, A0, NE, indexGas, indexCondensed, H0RT);

    % Equilibrium derivative respect pressure
    [dNi_p, dN_p] = equilibrium_dp(J, N, A0, NE, indexGas, indexCondensed);
end

function [dNi_T, dN_T] = equilibrium_dT(J, N, A0, NE, indexGas, indexCondensed, H0RT)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to temperature from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
    %     H0RT (float): Dimensionless enthalpy
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_T (float): Thermodynamic derivative of the moles of the species respect to temperature
    %     * dN_T (float):  Thermodynamic derivative of the moles of the mixture respect to temperature
    %
    % Example:
    %     [dNi_T, dN_T] = equilibrium_dT(J, N, A0, NE, ind, indexGas, indexCondensed, H0RT)

    % Definitions
    opts.SYM = true; % Options linsolve method: real symmetric

    % Initialization
    dNi_T = zeros(length(N), 1);

    % Construction of vector b
    b = update_vector_b(A0, N, indexGas, indexCondensed, H0RT);

    % Solve of the linear system J*x = b
    x = linsolve(J, b, opts);

    % Extract solution
    dpii_T = x(1:NE);
    dNi_T(indexCondensed) = x(NE+1:end-1);
    dN_T = x(end);

    % Compute remainder dNi_T (gas)
    dNi_T(indexGas) = H0RT(indexGas) + A0(indexGas, :) * dpii_T + dN_T;

    % NESTED FUNCTION
    function b = update_vector_b(A0, N, indexGas, indexCondensed, H0RT)
        % Compute vector b
        b = -[(sum(A0(indexGas, :) .* N(indexGas) .* H0RT(indexGas)))';...
              H0RT(indexCondensed); sum(N(indexGas) .* H0RT(indexGas))];
    end

end

function [dNi_p, dN_p] = equilibrium_dp(J, N, A0, NE, indexGas, indexCondensed)
    % Obtain thermodynamic derivative of the moles of the species and of the moles of the mixture
    % respect to pressure from a given composition [moles] at equilibrium
    %
    % Args:
    %     J (float): Matrix J to solve the linear system J*x = b
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     NE (float): Temporal total number of elements
    %     indexGas (float): Temporal index of gaseous species in the final mixture
    %     indexCondensed (float): Temporal index of condensed species in the final mixture
    %
    % Returns:
    %     Tuple containing
    %
    %     * dNi_p (float): Thermodynamic derivative of the moles of the species respect to pressure
    %     * dN_p (float):  Thermodynamic derivative of the moles of the mixture respect to pressure
    %
    % Example:
    %     [dNi_p, dN_p] = equilibrium_dp(J, N, A0, NE, ind, ind_nswt, ind_swt, ind_elem)

    % Definitions
    opts.SYM = true; % Options linsolve method: real symmetric

    % Initialization
    dNi_p = zeros(length(N), 1);

    % Construction of vector b
    b = update_vector_b(J, N, indexGas);

    % Solve of the linear system J*x = b
    x = linsolve(J, b, opts);

    % Extract solution
    dpii_p = x(1:NE);
    dNi_p(indexCondensed) = x(NE+1:end-1);
    dN_p = x(end);
    
    % Compute remainder dNi_T (gas)
    dNi_p(indexGas) = -1 + A0(indexGas, :) * dpii_p + dN_p;

    % NESTED FUNCTION
    function b = update_vector_b(J, N, ind_nswt)
        % Compute vector b
        b = J(:, end);
        b(end) = sum(N(ind_nswt));
    end

end
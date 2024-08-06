function [N, STOP, FLAG_ION] = equilibriumCheckIons(obj, N, A0, ind_E, indexGas, indexIons)
    % Check convergence of ionized species in the mixture
    %
    % Args:
    %     obj (EquilibriumSolver): EquilibriumSolver object
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     ind_E (float): Index of electron element
    %     indexGas (float): List of chemical species indices in gaseous phase
    %     indexIons (float): List of ionized chemical species indices 
    %
    % Returns:
    %     Tuple containing
    %
    %     * N (float): Mixture composition [mol]
    %     * STOP (float): Relative error in charge balance [-]
    %     * FLAG_ION (bool): Flag indicating if ionized species are present in the mixture
    
    % Initialization
    STOP = 0;

    % Check if there are ionized species
    if ~any(indexIons)
        FLAG_ION = false;
        return
    end
    
    % Update FLAG_ION
    FLAG_ION = true;

    % Get error in the electro-neutrality of the mixture
    [delta_ions, ~] = ionsFactor(N, A0, ind_E, indexGas, indexIons);
    
    % Reestimate composition of ionized species
    if abs(delta_ions) > obj.tolMultiplierIons
        [N, STOP] = recomputeIons(N, A0, ind_E, indexGas, indexIons, delta_ions, obj.tolMoles, obj.tolMultiplierIons, obj.itMaxIons);
    end
    
end

% SUB-PASS FUNCTIONS
function [N, STOP] = recomputeIons(N, A0, ind_E, indexGas, indexIons, delta_ions, TOL, TOL_pi, itMax)
    % Reestimate composition of ionized species
    %
    % Args:
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     ind_E (float): Index of electron element
    %     indexGas (float): List of gaseous species indices
    %     indexIons (float): List of ionized species indices
    %     delta_ions (float): Initial error in charge balance
    %     TOL (float): Tolerance for molar fraction
    %     TOL_pi (float): Tolerance for charge balance error
    %     itMax (int): Maximum number of iterations
    %
    % Returns:
    %     Tuple containing
    %
    %     * N (float): Updated mixture composition [mol]
    %     * STOP (float): Final relative error in charge balance
    
    % Initialization
    A0_ions = A0(indexIons, ind_E);
    STOP = 1;
    it = 0;

    % Reestimate composition of ionized species
    while STOP > TOL_pi && it < itMax
        % Update iteration
        it = it + 1;
        % Apply correction
        N(indexIons) = N(indexIons) .* exp(A0_ions * delta_ions);
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [delta_ions, ~] = ionsFactor(N, A0, ind_E, indexGas, indexIons);
        STOP = abs(delta_ions);
    end   
    
    Xi_ions = N(indexIons) / sum(N);

    % Set error to zero if molar fraction of ionized species are below tolerance
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end

end

function [delta, deltaN3] = ionsFactor(N, A0, ind_E, indexGas, indexIons)
    % Compute relaxation factor for ionized species
    %
    % Args:
    %     N (float): Mixture composition [mol]
    %     A0 (float): Stoichiometric matrix
    %     ind_E (float): Index of electron element
    %     indexGas (float): List of gaseous species indices
    %     indexIons (float): List of ionized species indices
    %
    % Returns:
    %     Tuple containing
    %
    %     * delta (float): Relaxation factor for charge balance
    %     * deltaN3 (float): Absolute sum of ion contributions
    
    % Return zero if there are no ionized species
    if ~any(indexIons)
        delta = [];
        deltaN3 = 0;
        return
    end
    
    % Compute relaxation factor for charge balance
    delta = -sum(A0(indexGas, ind_E) .* N(indexGas))/ ...
             sum(A0(indexGas, ind_E).^2 .* N(indexGas));
    deltaN3 = abs(sum(N(indexGas) .* A0(indexGas, ind_E)));
end
function [N, STOP, FLAG_ION] = equilibriumCheckIons(obj, N, A0, ind_E, indexGas, indexIons)
    % Check convergence of ionized species
    %
    % Args:
    %     obj (EquilibriumSolver): Equilibrium solver object
    %     N (float): Composition matrix [n_i, FLAG_CONDENSED_i] for the given temperature [K] and pressure [bar] at equilibrium
    %     A0 (float): Pressure [bar]
    %     ind_E (float): Index electron
    %     indexGas (flaot): List of chemical species indices in gaseous phase
    %     indexIons (float): List of ionized chemical species indices 
    %
    % Returns:
    %     Tuple containing
    %
    %     * N (float): Composition matrix [n_i, FLAG_CONDENSED_i] for the given temperature [K] and pressure [bar] at equilibrium
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
    
    % Initialization
    A0_ions = A0(indexIons, ind_E);
    STOP = 1;
    it = 0;

    % Reestimate composition of ionized species
    while STOP > TOL_pi && it < itMax
        % Update iteration
        it = it + 1;
        % Apply correction
        N(indexIons, 1) = N(indexIons, 1) .* exp(A0_ions * delta_ions);
        % Compute correction of the Lagrangian multiplier for ions divided by RT
        [delta_ions, ~] = ionsFactor(N, A0, ind_E, indexGas, indexIons);
        STOP = abs(delta_ions);
    end   
    
    Xi_ions = N(indexIons, 1) / sum(N(:, 1));

    % Set error to zero if molar fraction of ionized species are below tolerance
    if ~any(Xi_ions > TOL)
        STOP = 0;
    end

end

function [delta, deltaN3] = ionsFactor(N, A0, ind_E, indexGas, indexIons)
    % Compute relaxation factor for ionized species
    
    if ~any(indexIons)
        delta = [];
        deltaN3 = 0;
        return
    end

    delta = -sum(A0(indexGas, ind_E) .* N(indexGas, 1))/ ...
             sum(A0(indexGas, ind_E).^2 .* N(indexGas, 1));
    deltaN3 = abs(sum(N(indexGas, 1) .* A0(indexGas, ind_E)));
end
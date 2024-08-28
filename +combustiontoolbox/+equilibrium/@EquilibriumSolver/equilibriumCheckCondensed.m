function [indexCondensed, FLAG_CONDENSED, dL_dnj] = equilibriumCheckCondensed(A0, pi_i, W, indexCondensed, muRT, NC_max, FLAG_ONE, FLAG_RULE)
    % Check condensed species that may appear at chemical equilibrium
    %
    % Args:
    %     A0 (float): Stoichiometric matrix
    %     pi_i (float): Dimensionless Lagrange multiplier vector
    %     W (float): Molecular mass [kg/mol] vector
    %     indexCondensed (float): Index condensed species to be considered
    %     muRT (float): List of chemical species indices in gaseous phase
    %     NC_max (float): Maximum number of condensed species (Gibbs phase rule)
    %     FLAG_ONE (bool): Flag indicating to include condensed species in the system one by one
    %     FLAG_RULE (bool): Flag indicating to include condensed species in the system up to the maximum number of condensed species that satisfy the Gibbs phase rule
    %
    % Returns:
    %     Tuple containing
    %
    %     * indexCondensed (float): Index condensed species that may appear at chemical equilibrium
    %     * FLAG_CONDENSED (bool): Flag indicating additional condensed species that may appear at chemical equilibrium
    %     * dL_dnj (float): Vapor pressure test vector of the species that may appear at chemical equilibrium
    
    % Initialization
    FLAG_CONDENSED = false;

    % Checks
    if isempty(indexCondensed)
        dL_dnj = [];
        return
    end
    
    % Get length condensed species
    NC = length(indexCondensed);

    for i = NC:-1:1
        % Only check if there were atoms of the species in the initial
        % mixture
        if ~sum(A0(indexCondensed(i), :))
            continue
        end

        % Calculate dLdnj of the condensed species
        dL_dnj(i) = (muRT(indexCondensed(i)) - dot(pi_i, A0(indexCondensed(i), :))) / W(i);
    end
    
    % Get condensed species that may appear at chemical equilibrium
    FLAG = dL_dnj < 0;

    % Check if any condensed species have to be considered
    if ~sum(FLAG)
        indexCondensed = [];
        return
    end

    % Get index of the condensed species to be added to the system
    indexCondensed = indexCondensed(FLAG);
    dL_dnj = dL_dnj(FLAG);
    
    % Testing
    if FLAG_RULE
        [~, temp] = sort(dL_dnj);
        indexCondensed = indexCondensed(temp(1:NC_max));
    elseif FLAG_ONE
        [~, temp] = min(dL_dnj);
        indexCondensed = indexCondensed(temp);
    end

    % Update flag
    FLAG_CONDENSED = true;
end
function [indexCondensed, FLAG_CONDENSED, dL_dnj] = equilibriumCheckCondensed(A0, pi_i, W, indexCondensed, muRT, NC_max, FLAG_ONE, FLAG_RULE)
    % Check condensed species
    
    % Initialization
    FLAG_CONDENSED = false;

    % Get length condensed species
    NC = length(indexCondensed);

    for i = NC:-1:1
        % Only check if there were atoms of the species in the initial
        % mixture
        if ~sum(A0(indexCondensed(i), :))
            continue
        end

        % Calculate dLdnj of the condensed species
        dL_dnj(i) = (muRT(indexCondensed(i)) - dot(pi_i, A0(indexCondensed(i), :))) / W(i); % / W(i);
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
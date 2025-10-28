function h0 = getEnthalpyArray(listSpecies, T, DB)
    % Function that computes the vector of enthalpy for the given
    % set of species [J/mol]
    %
    % Args:
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database of `Species` objects with custom thermodynamic polynomials functions
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [J/mol]
    %
    % Example:
    %     h0 = getEnthalpyArray({'H2O', 'CO2'}, 298.15, DB)

    persistent cachedSpecies cachedH0curves cachedMap

    % Initialization
    if isempty(cachedH0curves)
        cachedSpecies = {};
        cachedH0curves = {};
        % Use containers.Map for O(1) lookup instead of O(n) find
        cachedMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    end

    % Definitions
    numSpecies = length(listSpecies);

    % Preallocate
    h0 = zeros(numSpecies, 1);

    % Get thermodynamic array
    for i = 1:numSpecies
        % Definitions
        species = listSpecies{i};

        % Check if species data is already cached using Map for O(1) lookup
        if isKey(cachedMap, species)
            index = cachedMap(species);
        else
            % Load species data and cache it
            cachedSpecies{end+1} = species;
            cachedH0curves{end+1} = DB.(species).h0curve;
            index = length(cachedSpecies);
            cachedMap(species) = index;
        end

        % Compute enthalpy [J/mol]
        h0(i) = cachedH0curves{index}(T);
    end

end

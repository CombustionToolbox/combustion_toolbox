function g0 = getGibbsEnergyArray(listSpecies, T, DB)
    % Function that computes the vector of Gibbs free energy for the given
    % set of species [J/mol]
    %
    % Args:
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     DB (struct): Database of `Species` objects with custom thermodynamic polynomials functions
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     g0 = getGibbsEnergyArray({'H2O', 'CO2'}, 298.15, DB)

    persistent cachedSpecies cachedG0curves cachedMap

    % Initialization
    if isempty(cachedG0curves)
        cachedSpecies = {};
        cachedG0curves = {};
        % Use containers.Map for O(1) lookup instead of O(n) find
        cachedMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    end

    % Definitions
    numSpecies = length(listSpecies);

    % Preallocate
    g0 = zeros(numSpecies, 1);

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
            cachedG0curves{end+1} = DB.(species).g0curve;
            index = length(cachedSpecies);
            cachedMap(species) = index;
        end

        % Compute Gibbs free energy [J/mol]
        g0(i) = cachedG0curves{index}(T);
    end

end

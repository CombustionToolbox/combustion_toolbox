function obj = updatePropertiesMatrixThermo(obj, listSpecies, T, index, FLAG_H0)
    % Function that computes the vector of enthalpy for the given
    % set of species [J/mol]
    %
    % Args:
    %     obj (ChemicalSystem): ChemicalSystem object
    %     listSpecies (cell): List of species
    %     T (float): Temperature [K]
    %     index (float): Vector with the indexes of the species to update the properties matrix
    %     FLAG_H0 (bool): Boolean to indicate if the enthalpy must be computed
    %
    % Returns:
    %     obj (ChemicalSystem): ChemicalSystem object with updated properties matrix
    %
    % Example:
    %     obj = updatePropertiesMatrixThermo(obj, {'H2O', 'CO2'}, 298.15, index, true)

    persistent cachedSpecies cachedCPcurves cachedH0curves cachedS0curves

    % Definitions
    numSpecies = length(listSpecies);

    % Initialization
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedCPcurves = {};
        cachedH0curves = {};
        cachedS0curves = {};
    end

    % Get thermodynamic array
    for i = 1:numSpecies
        % Definitions
        species = listSpecies{i};

        % Check if species data is already cached
        indexCache = find(strcmp(cachedSpecies, species), 1);

        % Load species data and cache it
        if isempty(indexCache)
            cachedSpecies{end+1} = species;
            cachedCPcurves{end+1} = obj.species.(species).cpcurve;
            cachedS0curves{end+1} = obj.species.(species).s0curve;
            cachedH0curves{end+1} = obj.species.(species).h0curve;

            indexCache = length(cachedSpecies);
        end

        % Compute thermodynamic properties
        if FLAG_H0
            obj.propertiesMatrix(index(i), obj.ind_hi) = cachedH0curves{indexCache}(T); % [J/mol]
        end

        obj.propertiesMatrix(index(i), obj.ind_cpi) = cachedCPcurves{indexCache}(T); % [J/(mol-K)]
        obj.propertiesMatrix(index(i), obj.ind_si) = cachedS0curves{indexCache}(T); % [J/(mol-K)]
    end

end
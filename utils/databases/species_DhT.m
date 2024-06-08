function DhT = species_DhT(species, T, DB)
    % Compute thermal enthalpy [J/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [J/mol]
    %
    % Example:
    %     DhT = species_DhT('H2O', 300, DB)

    persistent cachedSpecies;
    persistent cachedDHTcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedDHTcurves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, species), 1);
    if isempty(index)
        % Load species data and cache it
        DhTcurve = DB.(species).DhTcurve;
        cachedSpecies{end+1} = species;
        cachedDHTcurves{end+1} = DhTcurve;
    else
        % Retrieve cached data
        DhTcurve = cachedDHTcurves{index};
    end
    
    % Compute thermal enthalpy [J/mol]
    DhT = DhTcurve(T);
end

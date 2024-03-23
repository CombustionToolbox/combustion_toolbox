function DhT = species_DhT(species, T, DB)
    % Compute thermal enthalpy [kJ/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [kJ/mol]
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
    idx = find(strcmp(cachedSpecies, species), 1);
    if isempty(idx)
        % Load species data and cache it
        DhTcurve = DB.(species).DhTcurve;
        cachedSpecies{end+1} = species;
        cachedDHTcurves{end+1} = DhTcurve;
    else
        % Retrieve cached data
        DhTcurve = cachedDHTcurves{idx};
    end
    
    % Compute thermal enthalpy [kJ/mol]
    DhT = DhTcurve(T) / 1000;
end
function cP = species_cP(species, T, DB)
    % Compute specific heat at constant pressure [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     species (char): Chemical species
    %     T (float): Temperature [K]
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Returns:
    %     cP (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %
    % Example:
    %     cP = species_cP('H2O', 300, DB)

    persistent cachedSpecies;
    persistent cachedCPcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedCPcurves = {};
    end
    
    % Check if species data is already cached
    idx = find(strcmp(cachedSpecies, species), 1);
    if isempty(idx)
        % Load species data and cache it
        cPcurve = DB.(species).cPcurve;
        cachedSpecies{end+1} = species;
        cachedCPcurves{end+1} = cPcurve;
    else
        % Retrieve cached data
        cPcurve = cachedCPcurves{idx};
    end
    
    % Compute specific heat at constant pressure [J/(mol-K)]
    cP = cPcurve(T);
end

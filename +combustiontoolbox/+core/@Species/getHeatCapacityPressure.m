function cp = getHeatCapacityPressure(obj, T)
    % Compute specific heat at constant pressure [J/(mol-K)] of the species
    % at the given temperature [K] using piecewise cubic Hermite
    % interpolating polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     cp (float): Specific heat at constant pressure in molar basis [J/(mol-K)]
    %
    % Example:
    %     cp = getHeatCapacityPressure(obj, 300)
    
    persistent cachedSpecies;
    persistent cachedCPcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedCPcurves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        cpcurve = obj.cpcurve;
        cachedSpecies{end+1} = obj.name;
        cachedCPcurves{end+1} = cpcurve;
    else
        % Retrieve cached data
        cpcurve = cachedCPcurves{index};
    end
    
    % Compute specific heat at constant pressure [J/(mol-K)]
    cp = cpcurve(T);
end
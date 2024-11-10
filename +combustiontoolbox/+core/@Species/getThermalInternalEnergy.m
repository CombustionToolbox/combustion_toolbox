function DeT = getThermalInternalEnergy(obj, T)
    % Compute thermal internal energy [J/mol] of the species at the given
    % temperature [K] using piecewise cubic Hermite interpolating
    % polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     DeT (float): Thermal internal energy in molar basis [J/mol]
    %
    % Example:
    %     DeT = getThermalInternalEnergy(obj, 300)

    persistent cachedSpecies;
    persistent cachedDETcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedDETcurves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        DeTcurve = obj.DeTcurve;
        cachedSpecies{end+1} = obj.name;
        cachedDETcurves{end+1} = DeTcurve;
    else
        % Retrieve cached data
        DeTcurve = cachedDETcurves{index};
    end
    
    % Compute thermal internal energy [J/mol]
    DeT = DeTcurve(T);
end
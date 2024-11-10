function DhT = getThermalEnthalpy(obj, T)
    % Compute thermal enthalpy [J/mol] of the species at the given
    % temperature [K] using piecewise cubic Hermite interpolating
    % polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     DhT (float): Thermal enthalpy in molar basis [J/mol]
    %
    % Example:
    %     DhT = getThermalEnthalpy(obj, 300)

    persistent cachedSpecies;
    persistent cachedDHTcurves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedDHTcurves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        DhTcurve = obj.DhTcurve;
        cachedSpecies{end+1} = obj.name;
        cachedDHTcurves{end+1} = DhTcurve;
    else
        % Retrieve cached data
        DhTcurve = cachedDHTcurves{index};
    end
    
    % Compute thermal enthalpy [J/mol]
    DhT = DhTcurve(T);
end

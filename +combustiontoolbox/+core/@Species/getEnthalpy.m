function h0 = getEnthalpy(obj, T)
    % Compute enthalpy [J/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     h0 (float): Enthalpy in molar basis [J/mol]
    %
    % Example:
    %     h0 = getEnthalpy(obj, 300)
    
    persistent cachedSpecies;
    persistent cachedH0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedH0curves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        h0curve = obj.h0curve;
        cachedSpecies{end+1} = obj.name;
        cachedH0curves{end+1} = h0curve;
    else
        % Retrieve cached data
        h0curve = cachedH0curves{index};
    end
    
    % Compute enthalpy [J/mol]
    h0 = h0curve(T);
end
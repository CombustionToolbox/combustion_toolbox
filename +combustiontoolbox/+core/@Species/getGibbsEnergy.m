function g0 = getGibbsEnergy(obj, T)
    % Compute Gibbs energy [J/mol] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     g0 (float): Gibbs energy in molar basis [J/mol]
    %
    % Example:
    %     g0 = getGibbsEnergy(obj, 300)
    
    persistent cachedSpecies;
    persistent cachedG0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedG0curves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        g0curve = obj.g0curve;
        cachedSpecies{end+1} = obj.name;
        cachedG0curves{end+1} = g0curve;
    else
        % Retrieve cached data
        g0curve = cachedG0curves{index};
    end
    
    % Compute Gibbs energy [J/mol]
    g0 = g0curve(T);
end
function s0 = getEntropy(obj, T)
    % Compute entropy [J/(mol-K)] of the species at the given temperature [K]
    % using piecewise cubic Hermite interpolating polynomials and linear extrapolation
    %
    % Args:
    %     obj (Species): Species object
    %     T (float): Temperature [K]
    %
    % Returns:
    %     s0 (float): Entropy in molar basis [J/(mol-K)]
    %
    % Example:
    %     s0 = getEntropy(obj, 300)
    
    persistent cachedSpecies;
    persistent cachedS0curves;
    
    if isempty(cachedSpecies)
        cachedSpecies = {};
        cachedS0curves = {};
    end
    
    % Check if species data is already cached
    index = find(strcmp(cachedSpecies, obj.name), 1);
    if isempty(index)
        % Load species data and cache it
        s0curve = obj.s0curve;
        cachedSpecies{end+1} = obj.name;
        cachedS0curves{end+1} = s0curve;
    else
        % Retrieve cached data
        s0curve = cachedS0curves{index};
    end
    
    % Compute entropy [J/(mol-K)]
    s0 = s0curve(T);
end
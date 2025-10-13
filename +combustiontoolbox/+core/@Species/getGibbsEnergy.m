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
    
    persistent cachedID cachedG0curves
    
    % Get id
    id = obj.id_;

    % Preallocate cache if empty
    if isempty(cachedID)
        cachedID = zeros(0, 1, 'double');
        cachedG0curves = {};
    end
    
    % Check if species data is already cached
    index = find( cachedID == id, 1 );

    if isempty(index)
        % Load species data and cache it
        index = length(cachedID) + 1;
        cachedID(index) = id;
        cachedG0curves{index} = obj.g0curve;
    end

    % Compute Gibbs energy [J/mol]
    g0 = cachedG0curves{index}(T);
end
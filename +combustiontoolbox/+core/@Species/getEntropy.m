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
    
    persistent cachedID cachedS0curves
    
    % Get id
    id = obj.id_;

    % Preallocate cache if empty
    if isempty(cachedID)
        cachedID = zeros(0, 1, 'double');
        cachedS0curves = {};
    end
    
    % Check if species data is already cached
    index = find( cachedID == id, 1 );

    if isempty(index)
        % Load species data and cache it
        index = length(cachedID) + 1;
        cachedID(index) = id;
        cachedS0curves{index} = obj.s0curve;
    end

    % Compute entropy [J/(mol-K)]
    s0 = cachedS0curves{index}(T);
end
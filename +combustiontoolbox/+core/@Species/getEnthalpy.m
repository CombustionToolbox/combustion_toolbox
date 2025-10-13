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
    
    persistent cachedID cachedH0curves
    
    % Get id
    id = obj.id_;

    % Preallocate cache if empty
    if isempty(cachedID)
        cachedID = zeros(0, 1, 'double');
        cachedH0curves = {};
    end
    
    % Check if species data is already cached
    index = find( cachedID == id, 1 );

    if isempty(index)
        % Load species data and cache it
        index = length(cachedID) + 1;
        cachedID(index) = id;
        cachedH0curves{index} = obj.h0curve;
    end
    
    % Compute enthalpy [J/mol]
    h0 = cachedH0curves{index}(T);
end
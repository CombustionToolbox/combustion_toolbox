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
    
    persistent cachedID cachedCPcurves
    
    % Get id
    id = obj.id_;

    % Preallocate cache if empty
    if isempty(cachedID)
        cachedID = zeros(0, 1, 'double');
        cachedCPcurves = {};
    end
    
    % Check if species data is already cached
    index = find( cachedID == id, 1 );

    if isempty(index)
        % Load species data and cache it
        index = length(cachedID) + 1;
        cachedID(index) = id;
        cachedCPcurves{index} = obj.cpcurve;
    end

    % Compute specific heat at constant pressure [J/(mol-K)]
    cp = cachedCPcurves{index}(T);
end
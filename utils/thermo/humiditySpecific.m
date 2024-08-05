function value = humiditySpecific(T, p, humidityRelative)
    % Get the specific humidity of air [kg_w/kg_da] at a given temperature, pressure, and relative humidity
    %
    % Args:
    %     T (float): Temperature [K]
    %     p (float): Pressure [bar]
    %     humidityRelative (float): Relative humidity [%]
    %
    % Returns:
    %     value (float): Specific humidity of air [kg_w/kg_da]

    % Constants
    A = [-5.8002206e3, 1.3914993e0, -4.8640239e-2, 4.1764768e-5, -1.4452093e-8, 6.5459673e0];

    % Change pressure from [bar] to [Pa]
    p = p * 1e5;

    % Saturated vapor pressure [Pa]
    p_ws = exp( A(1)./T + A(2) + A(3) * T + A(4) * T.^2 + A(5) * T.^3 + A(6) * log(T) ); % Temperature in [K]

    % Vapor pressure [Pa]
    p_w = humidityRelative * p_ws * 1e-2;
    
    % Specific humidity [kg_w/kg_da]
    value = 0.622 * p_w ./ (p - p_w);
end

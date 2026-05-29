classdef EquationStateIdealGas < combustiontoolbox.core.EquationState
    % The :mat:func:`EquationStateIdealGas` class implements the ideal gas equation of state.
    %
    % The :mat:func:`EquationStateIdealGas` object can be used to compute the pressure,
    % molar volume, dimensionless volume derivatives, and thermodynamic departure functions
    % of a mixture assuming ideal gas behavior.
    %
    % Example:
    %      eos = EquationStateIdealGas();
    %
    % See also: :mat:func:`EquationState` and :mat:func:`Mixture`
    
    properties (Constant, Access = private)
        R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
    end

    methods (Access = public)

        function pressure = getPressure(obj, temperature, molarVolume, ~, ~, varargin)
            % Compute pressure [Pa] using the ideal gas law.
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %
            % Returns:
            %     pressure (float): Pressure of the mixture [Pa]
            %
            % Example:
            %     P = getPressure(obj, 300, 0.024)

            % Compute pressure [Pa]
            pressure = (obj.R0 .* temperature) ./ molarVolume;
        end
        
        function molarVolume = getVolume(obj, temperature, pressure, ~, ~, varargin)
            % Compute molar volume [m3/mol] using the ideal gas law.
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %
            % Returns:
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %
            % Example:
            %     V = getVolume(obj, 300, 1e5)

            % Compute molar volume [m3/mol]
            molarVolume = obj.R0 * temperature ./ pressure;
        end

        function [dPdV_T, dPdT_V] = getPressureDerivativesDimensional(obj, temperature, ~, molarVolume, ~, ~, varargin)
            % Compute dimensional partial pressure derivatives for an ideal gas
            %
            % Args:
            %     obj (EquationStateIdealGas): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %
            % Returns:
            %     Tuple containing
            %
            %     * dPdV_T (float): Partial derivative of pressure with respect to volume at constant temperature [Pa/(m3/mol)]
            %     * dPdT_V (float): Partial derivative of pressure with respect to temperature at constant volume [Pa/K]
            %
            % Example:
            %     [dPdV_T, dPdT_V] = getPressureDerivativesDimensional(obj, 300, 1e5, 0.024)

            % Compute dimensional pressure derivatives for Ideal Gas: P = R0*T/V
            dPdV_T = -(obj.R0 * temperature) / (molarVolume^2);
            dPdT_V = obj.R0 / molarVolume;
        end

        function [dVdT_p, dVdp_T] = getVolumeDerivatives(~, ~, ~, ~, ~, ~, varargin)
            % Compute dimensionless volume derivatives for the mixture assuming frozen chemistry.
            %
            % Returns:
            %     Tuple containing
            %
            %     * dVdT_p (float): Logarithmic derivative of volume with respect to temperature at constant pressure [-]
            %     * dVdp_T (float): Logarithmic derivative of volume with respect to pressure at constant temperature [-]
            %
            % Example:
            %     [dVdT_p, dVdp_T] = getVolumeDerivatives(obj)

            % For an ideal gas (V = RT/p), the dimensionless derivatives are exactly 1 and -1.
            dVdT_p = 1;
            dVdp_T = -1;
        end

        function [heatCapacityPressureDeparture, enthalpyDeparture, entropyDeparture] = getDepartureFunctions(~, ~, ~, ~, ~, ~, varargin)
            % Compute thermodynamic departure functions for an ideal gas. For an ideal gas, all departure functions are zero.
            %
            % Returns:
            %     Tuple containing
            %
            %     * heatCapacityPressureDeparture (float): Heat capacity at constant pressure departure [J/(mol-K)]
            %     * enthalpyDeparture (float): Enthalpy departure [J/mol]
            %     * entropyDeparture (float): Entropy departure [J/(mol-K)]
            %
            % Example:
            %     [dcp, dh, ds] = getDepartureFunctions(obj)

            heatCapacityPressureDeparture = 0;
            enthalpyDeparture = 0;
            entropyDeparture = 0;
        end

        function temperature = getTemperature(obj, pressure, molarVolume, varargin)
            % Compute temperature [K] given the pressure and molar volume using the ideal gas law.
            %
            % Args:
            %     pressure (float): Target pressure of the mixture [Pa]
            %     molarVolume (float): Target molar volume of the mixture [m3/mol]
            %
            % Returns:
            %     temperature (float): Computed temperature of the mixture [K]

            temperature = (pressure .* molarVolume) ./ obj.R0;
        end

    end

end

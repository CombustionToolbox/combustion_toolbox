classdef (Abstract) EquationState < handle
    % The :mat:func:`EquationState` class is an abstract class that defines
    % the interface for computing the pressure and molar volume of a
    % mixture using a specified equation of state.
    %
    % Subclasses must implement the following abstract methods:
    %     * getPressure(obj,temperature, molarVolume, varargin)
    %     * getVolume(obj, temperature, pressure, varargin)
    %     * getDepartureFunctions(obj, temperature, pressure, varargin)
    %     * getVolumeDerivatives(obj, temperature, pressure, molarVolume, varargin)
    %
    % See also: :mat:func:`EquationStateIdealGas`, :mat:func:`Mixture`

    methods (Abstract)
        % Compute pressure [Pa] given the temperature and molar volume.
        pressure = getPressure(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
        
        % Compute molar volume [m3/mol] given the temperature and pressure.
        molarVolume = getVolume(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)

        % Compute thermodynamic departure functions
        [heatCapacityPressureDeparture, enthalpyDeparture, entropyDeparture] = getDepartureFunctions(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)

        % Compute dimensionless volume derivatives: (dlnV/dlnT)_p and (dlnV/dlnP)_T
        [dlnVdT_p, dlnVdp_T] = getVolumeDerivatives(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
    end

    methods (Access = public, Static)

        function Z = getCompressibilityFactor(temperature, pressure, molarVolume)
            % Compute compressibility factor
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %
            % Returns:
            %     Z (float): Compressibility factor [-]
            %
            % Example:
            %     Z = EquationState.getCompressibilityFactor(300, 1e5, 0.024)

            % Definitions
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            
            % Compute compressibility factor
            Z = (pressure * molarVolume) / (R0 * temperature);
        end

    end
    
end
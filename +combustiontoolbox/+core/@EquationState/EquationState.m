classdef (Abstract) EquationState < handle
    % The :mat:func:`EquationState` class is an abstract class that defines
    % the interface for computing the pressure and molar volume of a
    % mixture using a specified equation of state.
    %
    % Subclasses must implement the following abstract methods:
    %     * getPressure(obj,temperature, molarVolume, molarFractions, chemicalSystem, varargin)
    %     * getVolume(obj, temperature, pressure, molarFractions, chemicalSystem, varargin)
    %     * getDepartureFunctions(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
    %     * getPressureDerivativesDimensional(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
    %
    % See also: :mat:func:`EquationStateIdealGas`, :mat:func:`Mixture`

    methods (Abstract)
        % Compute pressure [Pa] given the temperature and molar volume.
        pressure = getPressure(obj, temperature, molarVolume, molarFractions, chemicalSystem, varargin)
        
        % Compute molar volume [m3/mol] given the temperature and pressure.
        molarVolume = getVolume(obj, temperature, pressure, molarFractions, chemicalSystem, varargin)

        % Compute thermodynamic departure functions
        [heatCapacityPressureDeparture, enthalpyDeparture, entropyDeparture] = getDepartureFunctions(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)

        % Compute dimensional pressure derivatives: (dP/dV)_T and (dP/dT)_V
        [dPdV_T, dPdT_V] = getPressureDerivativesDimensional(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin);
    end

    methods (Access = public)

        function [dPdV_T, dPdT_V] = getPressureDerivatives(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute dimensionless (logarithmic) partial pressure derivatives for the mixture
            %
            % Args:
            %     obj (EquationState): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * dPdV_T (float): Logarithmic derivative of pressure with respect to volume at constant temperature [-]
            %     * dPdT_V (float): Logarithmic derivative of pressure with respect to temperature at constant volume [-]

            % Get dimensional pressure derivatives from the specific EoS implementation
            [dPdV_T, dPdT_V] = obj.getPressureDerivativesDimensional(temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin{:});
            
            % Convert to dimensionless (logarithmic) pressure derivatives
            % (dlnP/dlnV)_T = (V/P) * (dP/dV)_T
            % (dlnP/dlnT)_V = (T/P) * (dP/dT)_V
            dPdV_T = (molarVolume / pressure) * dPdV_T;
            dPdT_V = (temperature / pressure) * dPdT_V;
        end

        function [dVdT_p, dVdp_T] = getVolumeDerivatives(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute dimensionless (logarithmic) volume derivatives for the mixture assuming frozen chemistry
            %
            % Args:
            %     obj (EquationStatePengRobinson): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * dVdT_p (float): Logarithmic derivative of volume with respect to temperature at constant pressure [-]
            %     * dVdp_T (float): Logarithmic derivative of volume with respect to pressure at constant temperature [-]
            %
            % Example:
            %     [dVdT_p, dVdp_T] = getVolumeDerivatives(obj, 300, 1e5, 0.024, [0.5, 0.5], chemicalSystem)

            % Compute dimensional pressure derivatives
            [dPdV_T, dPdT_V] = obj.getPressureDerivativesDimensional(temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin{:});
            
            % Convert to volume derivatives using the chain rule:
            %
            %   (dlnV/dlnT)_p = -(T/V) * (dP/dT)_V / (dP/dV)_T 
            %   (dlnV/dlnP)_T =  (P/V) * (1 / (dP/dV)_T)
            dVdT_p = (temperature / molarVolume) * (-dPdT_V / dPdV_T);
            dVdp_T = (pressure / molarVolume) * (1 / dPdV_T);
        end

        function [dVdT_p, dVdp_T] = getVolumeDerivativesDimensional(obj, temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin)
            % Compute dimensional volume derivatives for the mixture assuming frozen chemistry
            %
            % Args:
            %     obj (EquationState): Equation of state object
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Returns:
            %     Tuple containing
            %
            %     * dVdT_p (float): Partial derivative of volume with respect to temperature at constant pressure [m3/(mol-K)]
            %     * dVdp_T (float): Partial derivative of volume with respect to pressure at constant temperature [m3/(mol-Pa)]
            %
            % Example:
            %     [dVdT_p, dVdp_T] = getVolumeDerivativesDimensional(obj, 300, 1e5, 0.024, [0.5, 0.5], chemicalSystem)

            % Get dimensional pressure derivatives from the specific EoS implementation
            [dPdV_T, dPdT_V] = obj.getPressureDerivativesDimensional(temperature, pressure, molarVolume, molarFractions, chemicalSystem, varargin{:});
            
            % Compute dimensional volume derivatives
            % (dV/dT)_p = -(dP/dT)_V / (dP/dV)_T
            % (dV/dp)_T = 1 / (dP/dV)_T
            dVdT_p = -dPdT_V / dPdV_T;
            dVdp_T = 1 / dPdV_T;
        end

        function temperature = getTemperature(obj, pressure, molarVolume, molarFractions, chemicalSystem, temperatureGuess, varargin)
            % Compute temperature [K] given the pressure and molar volume using a numerical root-finder
            %
            % Args:
            %     obj (EquationState): Equation of state object
            %     pressure (float): Pressure of the mixture [Pa]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
            %     molarFractions (float): Molar fractions of the species in the mixture
            %     chemicalSystem (ChemicalSystem): Chemical system object containing species data
            %
            % Optional Args:
            %     temperatureGuess (float): Initial guess for the temperature [K] (default: 300 K)
            %
            % Returns:
            %     temperature (float): Temperature of the mixture [K]

            % Definitions
            DefaultTemperatureGuess = 300; % [K]

            if nargin < 6 || isempty(temperatureGuess)
                temperatureGuess = DefaultTemperatureGuess;
            end

            % Set options for fzero
            options = optimset('Display', 'off');
            
            % Determine output size and expand scalars if necessary
            if isscalar(pressure)
                sizeStates = size(molarVolume);
                pressure = repmat(pressure, sizeStates);
            elseif isscalar(molarVolume)
                sizeStates = size(pressure);
                molarVolume = repmat(molarVolume, sizeStates);
            else
                sizeStates = size(pressure);
                assert(isequal(size(pressure), size(molarVolume)), 'pressure and molarVolume must have the same size unless one is scalar.');
            end

            numStates = prod(sizeStates);

            % Vectorize for faster linear indexing
            pressure = pressure(:);
            molarVolume = molarVolume(:);
            
            % Loop over all target states
            for i = numStates:-1:1
                pTarget = pressure(i);
                vTarget = molarVolume(i);

                % Objective function: difference between EOS pressure and target pressure
                pressure_error = @(T) obj.getPressure(T, vTarget, molarFractions, chemicalSystem, varargin{:}) - pTarget;
                
                % Solve for T using fzero
                temperature(i) = fzero(pressure_error, temperatureGuess, options);

                % Update guess for next iteration
                temperatureGuess = temperature(i);
            end

            % Reshape output to match input size
            temperature = reshape(temperature, sizeStates);
        end

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
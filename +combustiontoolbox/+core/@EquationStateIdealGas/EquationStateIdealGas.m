classdef EquationStateIdealGas < combustiontoolbox.core.EquationState
    % The :mat:func:`EquationStateIdealGas` class implements the ideal gas law.
    %
    % The :mat:func:`EquationStateIdealGas` object can be used to compute
    % the pressure and molar volume of a mixture assuming ideal gas behavior.
    %
    % Example:
    %      eos = EquationStateIdealGas();
    %
    % See also: :mat:func:`EquationState` and :mat:func:`Mixture`
    
    methods (Access = public, Static)

        function pressure = getPressure(temperature, molarVolume, varargin)
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

            % Definitions
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]

            % Compute pressure [Pa]
            pressure = (R0 .* temperature) ./ molarVolume;
        end
        
        function molarVolume = getVolume(temperature, pressure, varargin)
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

            % Definitions
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]

            % Compute molar volume [m3/mol]
            molarVolume = R0 * temperature ./ pressure;
        end

    end

end

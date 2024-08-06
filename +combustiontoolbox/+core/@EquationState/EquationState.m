classdef EquationState < handle
    % The :mat:func:`EquationState` class is used to compute the pressure and molar volume of a mixture using the equation of state specified by the user.
    %
    % The :mat:func:`EquationState` object can be initialized as follows:
    %
    %       equationState = EquationState('idealgas')
    %
    % This creates an instance of the :mat:func:`EquationState` class and initializes it with the
    % equation of state specified by the user.
    %
    % See also: :mat:func:`Mixture`
    
    properties
        equationOfState
        pressureFunction
        volumeFunction
    end
    
    methods

        function obj = EquationState(equationOfState)
            % Constructor
            obj.equationOfState = equationOfState;
            obj = obj.setEquationState();
        end
        
        function obj = setEquationState(obj)
            % Compute pressure using the equation of state specified by the user
            switch lower(obj.equationOfState)
                case 'idealgas'
                    obj.pressureFunction = @obj.getPressureIdeal;
                    obj.volumeFunction = @obj.getVolumeIdeal;
                otherwise
                    error('Invalid equation of state');
            end
        end

        function pressure = getPressure(obj, temperature, molarVolume, varargin)
            % Compute pressure [Pa] using the defined equation of state
            pressure = obj.pressureFunction(temperature, molarVolume, varargin{:});
        end

        function volume = getVolume(obj, temperature, pressure, varargin)
            % Compute molar volume [m3/mol_gas] using the defined equation of state
            volume = obj.volumeFunction(temperature, pressure, varargin{:});
        end

    end

    methods (Access = private, Static)
        
        function pressure = getPressureIdeal(temperature, molarVolume, varargin)
            % Compute pressure considering ideal Equation of State (EoS)
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     molarVolume (float): Molar volume of the mixture [m3/mol_gas]
            % 
            % Returns:
            %     pressure (float): Pressure of the mixture [Pa]
            
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
        
            pressure = (R0 .* temperature) ./ molarVolume;
        end

        function molarVolume = getVolumeIdeal(temperature, pressure, varargin)
            % Compute molar volume considering ideal Equation of State (EoS)
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     pressure (float): Pressure of the mixture [Pa]
            % 
            % Returns:
            %     molarVolume (float): Molar volume of the mixture [m3/mol_gas]
            
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            
            molarVolume = R0 * temperature ./ pressure;
        end

    end

end
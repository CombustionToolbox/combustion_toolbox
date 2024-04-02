classdef EquationState < handle
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
            % Compute molar volume [m3/mol] using the defined equation of state
            volume = obj.volumeFunction(temperature, pressure, varargin{:});
        end

    end

    methods (Access = private, Static)
        
        function pressure = getPressureIdeal(temperature, molarVolume, varargin)
            % Compute pressure considering ideal Equation of State (EoS)
            %
            % Args:
            %     temperature (float): Temperature of the mixture [K]
            %     molarVolume (float): Molar volume of the mixture [m3/mol]
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
            %     molarVolume (float): Pressure of the mixture [Pa]
            % 
            % Returns:
            %     V (float): Molar volume of the mixture [m3/mol]
            
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            
            molarVolume = R0 * temperature ./ pressure;
        end

    end

end
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

        function pressure = getPressure(obj, molesGas, temperature, volume, varargin)
            % Compute pressure [Pa] using the equation of state specified by the user            
            pressure = obj.pressureFunction(molesGas, temperature, volume, varargin{:});
        end

        function volume = getVolume(obj, molesGas, temperature, pressure, varargin)
            % Compute molar volume [m3/mol] using the equation of state specified by the user
            volume = obj.volumeFunction(molesGas, temperature, pressure, varargin{:});
        end

    end

    methods (Access = private, Static)
        
        function p = getPressureIdeal(n, T, v, varargin)
            % Compute pressure considering ideal Equation of State (EoS)
            %
            % Args:
            %     n (float): Number of moles of the mixture in gaseous phase [mol]
            %     T (float): Temperature of the mixture [K]
            %     v (float): Volume of the mixture [m3]
            % 
            % Returns:
            %     p (float): Pressure of the mixture [Pa]
            
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
        
            p = (n .* R0 .* T) ./ v;
        end

        function V = getVolumeIdeal(T, p, varargin)
            % Compute molar volume considering ideal Equation of State (EoS)
            %
            % Args:
            %     T (float): Temperature of the mixture [K]
            %     p (float): Pressure of the mixture [Pa]
            % 
            % Returns:
            %     V (float): Molar volume of the mixture [m3/mol]
            
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            
            V = R0 * T ./ p;
        end

    end

end
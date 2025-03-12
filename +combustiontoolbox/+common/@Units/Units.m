classdef Units < handle
    % Class with conversion factors between different units
    
    properties (Constant)
        % Pressure conversion factors
        atm2bar = @(x) x * 1.01325     % Atmospheres to bar
        bar2atm = @(x) x * 1.01325^-1  % Bar to atmospheres
        atm2Pa  = @(x) x * 101325      % Atmospheres to Pascales
        Pa2atm  = @(x) x * 101325^-1   % Pascales to atmospheres
        bar2Pa  = @(x) x * 1e5         % Bar to Pascals
        Pa2bar  = @(x) x * 1e-5        % Pascals to bar
        % Temperature conversion factors
        K2C     = @(x) x - 273.15      % Kelvin to degrees Celsius
        C2K     = @(x) x + 273.15      % Degrees Celsius to Kelvin
        F2C     = @(x) (x - 32) * 5/9  % Fahrenheit to degrees Celsius
        K2F     = @(x) (x - 273.15) * 9/5 + 32 % Kelvin to Fahrenheit
        % Mass conversion factors
        kg2lbs  = @(x) x * 2.20462     % Kilograms to pounds
        lbs2kg  = @(x) x * 0.453592    % Pounds to kilograms
        kg2g    = @(x) x * 1e3         % Kilograms to grams
        g2kg    = @(x) x * 1e-3        % Grams to kilograms
        % Volume conversion factors
        m32ft3  = @(x) x * 35.3147     % Cubic meters to cubic feet
        ft32m3  = @(x) x * 35.3147^-1  % Cubic feet to cubic meters
        m32L    = @(x) x * 1e3         % Cubic meters to liters
        L2m3    = @(x) x * 1e-3        % Liters to cubic meters
        ft32L   = @(x) x * 28.3168     % Cubic feet to liters
        L2ft3   = @(x) x * 28.3168^-1  % Liters to cubic feet
    end

    methods (Static)

        function value_out = convert(value_in, unit_in, unit_out)
            % Convert a value from one unit to another
            %
            % Args:
            %     value_in (double): Value to convert
            %     unit_in (char): Unit of the input value
            %     unit_out (char): Unit of the output value
            %
            % Returns:
            %     value_out: Value converted to the output unit
            %
            % Examples:
            %     * Units.convert(1, 'atm', 'bar')
            %     * Units.convert(273.15, 'K', 'C')
            %     * Units.convert(1, 'kg', 'lbs')
            %     * Units.convert(1, 'm3', 'ft3')
            %     * Units.convert(1, 'm3', 'L')
            
            % Import packages
            import combustiontoolbox.common.Units
            
            % Get the conversion key property name
            conversionKey = [unit_in, '2', unit_out];
            
            % Check conversion key exist
            assert(isprop(Units, conversionKey), 'Conversion from %s to %s is not defined.', unit_in, unit_out); 

            % Get the conversion key value
            conversion = Units.(conversionKey);
            
            % Convert the value
            value_out = conversion(value_in);
        end

        function moles = convertWeightPercentage2moles(listSpecies, weightPercentage, database)
            % Convert weight percentage (wt%) to moles
            %
            % Args:
            %     listSpecies (cell): List of species
            %     weightPercentage (float): Weight percentage of the species [%]
            %     database (Database): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %
            % Returns:
            %     moles (float): Number of moles [mol]
            %
            % Example:
            %     moles = Units.WeightPercentage2moles({'H2O', 'CO2'}, [50, 50], database)
        
            % Check if value is a cell
            if ~iscell(listSpecies)
                listSpecies = {listSpecies};
            end
        
            % Get molecular weight [g] of the species
            W = getProperty(database, listSpecies, 'W') * 1e3;
            
            % Convert weight percentage (wt%) to moles
            moles = weightPercentage ./ W;
        end

        function velocity = convertData2VelocityField(data)
            % Convert the input to a VelocityField object
            %
            % Args:
            %     data: Either a VelocityField object, a struct, or a 4D matrix
            %
            % Returns:
            %     velocity (VelocityField): VelocityField object
            
            % Import packages
            import combustiontoolbox.turbulence.VelocityField
            
            if isa(data, 'VelocityField')
                % Input is already a VelocityField object
                velocity = data;
                return
            end

            if isstruct(data) && all(isfield(data, {'u', 'v', 'w'}))
                % Input is a struct; convert to VelocityField
                velocity = VelocityField(data.u, data.v, data.w);
                return;
            end

            if ndims(data) == 4 && size(data, 4) == 3
                % Input is a 4D matrix; convert to VelocityField
                velocity = VelocityField(data(:, :, :, 1), ...
                                         data(:, :, :, 2), ...
                                         data(:, :, :, 3));
                return;
            end
            
            % Unsupported format
            error('Unsupported velocity input format. Expected a VelocityField object, a struct, or a 4D matrix.');
        end

    end
    
end
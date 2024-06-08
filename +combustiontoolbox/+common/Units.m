classdef Units < handle
    % Class with conversion factors between different units
    
    properties (Constant)
        % Pressure
        atm2bar = 1.01325      % Conversion factor from atm to bar
        bar2atm = 1.01325^-1   % Conversion factor from bar to atm
        atm2Pa  = 101325       % Conversion factor from atm to Pa
        Pa2atm  = 101325^-1    % Conversion factor from Pa to atm
        bar2Pa  = 1e5          % Conversion factor from bar to Pa
        Pa2bar  = 1e-5         % Conversion factor from Pa to bar
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
            % Example:
            %     Units.convert(1, 'atm', 'bar')
            
            % Import packages
            import combustiontoolbox.common.Units
            
            % Get the conversion factor property name
            conversion_factor_name = [unit_in,'2',unit_out];
            
            % Check conversion factor exist
            assert(isprop(Units, conversion_factor_name), 'Conversion from %s to %s is not defined.', unit_in, unit_out); 

            % Get the conversion factor
            conversion_factor = Units.(conversion_factor_name);
            
            % Convert the value
            value_out = value_in * conversion_factor;
        end

        function moles = convertWeightPercentage2moles(listSpecies, weight_percentage, DB)
            % Convert weight percentage (wt%) to moles
            %
            % Args:
            %     listSpecies (cell): List of species
            %     weight_percentage (float): Weight percentage of the species [%]
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %
            % Returns:
            %     moles (float): Number of moles [mol]
            %
            % Example:
            %     moles = Units.WeightPercentage2moles({'H2O', 'CO2'}, [50, 50], DB)
        
            % Check if value is a cell
            if ~iscell(listSpecies)
                listSpecies = {listSpecies};
            end
        
            % Get molecular weight of the species
            W = set_prop_DB(listSpecies, 'mm', DB);
            
            % Convert weight percentage (wt%) to moles
            moles = weight_percentage ./ W;
        end

    end
    
end
classdef Export < handle
    % Class to export and save data to different formats (Excel, CSV, and mat)

    properties
        format     % Format to save data (Excel, CSV, or mat)
        filename   % Filename to save data
        FLAG_PROMT % Flag to show promt to save data
        compositionUnit = 'kg'; % Composition unit (kg or mol)
        compositionVariable = 'X'; % Composition variable (molar fraction)
    end

    methods
        function obj = Export(varargin)
            % Class constructor

            % Default values
            defaultFormat = '.xls';
            defaultFilename = 'data';
            defaultFLAG_PROMT = true;

            % Create input parser
            ip = inputParser;
            addParameter(ip, 'format', defaultFormat, @ischar);
            addParameter(ip, 'filename', defaultFilename, @ischar);
            addParameter(ip, 'FLAG_PROMT', defaultFLAG_PROMT, @islogical);
            parse(ip, varargin{:});

            % Set properties
            obj.format = ip.Results.format;
            obj.filename = ip.Results.filename;
            obj.FLAG_PROMT = ip.Results.FLAG_PROMT;
        end

        function export(obj, mixArray, varargin)
            % Export results of reactants (mix1) and products (mix2) into a .xls
            % file
            %
            % Args:
            %     mixArray (Mixture): Array of mixture objects
            %
            % Optional Args:
            %     * mixArray_i (Mixture): Additional arrayy of mixture objects
            %
            % Examples:
            %     * exportResults(mixArray1)
            %     * exportResults(mixArray1, mixArray2)
            %     * exportResults(mixArray1, mixArray2, mixArray3)

            % Unpack additional inputs
            varargin = [{mixArray}, varargin];

            % Definitions
            numCases = length(varargin);
            equivalenceRatio = [mixArray.equivalenceRatio];
            sheetName = mixArray(1).problemType;

            % Format and save data
            for i = numCases:-1:1
                % Definitions
                mixArray = varargin{i};
                listSpecies = mixArray(1).chemicalSystem.listSpecies;

                % Format data
                data(:, :, i) = obj.generateMatrix(mixArray, listSpecies, equivalenceRatio);
            end

            % Show promt?
            if obj.FLAG_PROMT
                % Show dialog to save file
                [filename, filepath] = uiputfile(sprintf('*%s', obj.format));
            else
                filename = obj.filename;
                filepath = [];
            end

            % Export the figure if the user did not cancel
            if ~ischar(filename)
                return
            end

            % Assign filename to the path
            filename = [filepath filename];

            % Export data
            switch lower(obj.format)
                case {'.xls', 'xls', 'excel', 'spreadsheet'}

                    for i = numCases:-1:1
                        writecell(data(:, :, i), filename, 'Sheet', sprintf('Mixture%d-%s', i, sheetName));
                    end
                    
                case {'.mat', 'mat', 'matlab'}
                    save(filename, 'data');
            end

        end

    end

    methods (Access = public, Static)
        function data = generateMatrix(mixArray, listSpecies, equivalenceRatio)
            % Construct a cell with the thermodynamic data of the given mixture
            %
            % Args:
            %     mixArray (Mixture): Mixture object
            %     listSpecies (cell): List of species
            %     equivalenceRatio (float): Vector of equivalence ratio
            %
            % Returns:
            %     data (cell): Cell array with the thermodynamic data of the given mixture
            % 
            % Example:
            %     data(mixArray, {'CH4', 'O2', 'CO2', 'H2O'}, [0.5, 1.0, 1.5])
            
            % Definitions
            numCases = length(mixArray);
            numProperties = 19;
            numSpecies = numel(listSpecies);
            compositionUnit = obj.compositionUnit;
            compositionVariable = obj.compositionVariable;
            FLAG_VELOCITY = ~isempty(mixArray(1).u);
        
            % Get data phases species
            phase = mixArray(1).phase';
            
            % Set phase labels
            phaseLabels = cell(1, numSpecies);
            phaseLabels(phase > 0) = {'condensed'};
            phaseLabels(phase == 0) = {'gas'};
            
            % Headers - thermodynamic properties
            data{1, 1} = 'phi';      data{2, 1} = '[-]';
            data{1, 2} = 'T';        data{2, 2} = '[K]';
            data{1, 3} = 'p';        data{2, 3} = '[bar]';
            data{1, 4} = 'rho';      data{2, 4} = '[kg/m^3]';
            data{1, 5} = 'h';        data{2, 5} = ['[kJ/', compositionUnit,']'];
            data{1, 6} = 'e';        data{2, 6} = ['[kJ/', compositionUnit,']'];
            data{1, 7} = 's';        data{2, 7} = ['[kJ/(', compositionUnit,'·K)]'];
            data{1, 8} = 'g';        data{2, 8} = ['[kJ/', compositionUnit,']'];
            data{1, 9} = 'N';        data{2, 9} = '[mol]';
            data{1, 10} = 'W';       data{2, 10} = '[g/mol]';
            data{1, 11} = 'MW';      data{2, 11} = '[g/mol]';
            data{1, 12} = 'm';       data{2, 12} = '[kg]';
            data{1, 13} = 'dVdT_p';  data{2, 13} = '[-]';
            data{1, 14} = 'dVdp_T';  data{2, 14} = '[-]';
            data{1, 15} = 'cp';      data{2, 15} = ['[kJ/(', compositionUnit,'·K)]'];
            data{1, 16} = 'gamma';   data{2, 16} = '[-]';
            data{1, 17} = 'gamma_s'; data{2, 17} = '[-]';
            data{1, 18} = 'sound';   data{2, 18} = '[m/s]';
            data{1, 19} = 'u';       data{2, 19} = '[m/s]';
            
            % Headers - chemical species
            for j = numSpecies:-1:1
                data{1, numProperties + j} = sprintf('%s_%s (%s)', compositionVariable, listSpecies{j}, phaseLabels{j});
                data{2, numProperties + j} = '[-]';
            end
            
            % Assign data
            for i = numCases:-1:1
                % Equivalence ratio
                data{i + 2, 1} = equivalenceRatio(i);
                % Temperature
                data{i + 2, 2} = temperature(mixArray(i));
                % Pressure
                data{i + 2, 3} = pressure(mixArray(i));
                % Density
                data{i + 2, 4} = density(mixArray(i));
                % Specific enthalpy
                data{i + 2, 5} = enthalpy_mass(mixArray(i));
                % Specific internal energy
                data{i + 2, 6} = intEnergy_mass(mixArray(i));
                % Specific Gibbs energy
                data{i + 2, 7} = gibbs_mass(mixArray(i));
                % Specific entropy
                data{i + 2, 8} = entropy_mass(mixArray(i));
                % Total number of moles
                data{i + 2, 9} = mixArray(i).N;
                % Molecular weight (gases)
                data{i + 2, 10} = MolecularWeight(mixArray(i));
                % Mean molecular weight
                data{i + 2, 11} = meanMolecularWeight(mixArray(i));
                % Mass
                data{i + 2, 12} = mass(mixArray(i));
                % Thermodynamic derivative at constant pressure
                data{i + 2, 13} = mixArray(i).dVdT_p;
                % Thermodynamic derivative at constant temperature
                data{i + 2, 14} = mixArray(i).dVdp_T;
                % Specific heat at constant pressure
                data{i + 2, 15} = cp_mass(mixArray(i));
                % Specific heat ratio
                data{i + 2, 16} = adiabaticIndex(mixArray(i));
                % Adiabatic index
                data{i + 2, 17} = adiabaticIndex_sound(mixArray(i));
                % Sound velocity
                data{i + 2, 18} = soundspeed(mixArray(i));
                % Shock velocity
                if FLAG_VELOCITY
                    data{i + 2, 19} = velocity_relative(mixArray(i));
                end
            
                % Molar fractions
                for j = numSpecies:-1:1
                    data{i + 2, numProperties + j} = mixArray(i).Xi(j);
                end
                
            end
        
        end

    end

end
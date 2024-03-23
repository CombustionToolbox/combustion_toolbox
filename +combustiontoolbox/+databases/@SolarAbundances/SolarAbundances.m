classdef SolarAbundances
    % SolarAbundances class to read solar abundances and compute the initial molar composition
    % of the mixture 
    %
    % Optional Args:
    %     filename (char): Filename with the solar elemental abundances
    %
    % Examples:
    %     * DB = combustiontoolbox.databases.SolarAbundances;
    %     * moles = DB.abundances2moles({'H', 'He', 'C', 'N', 'O', 'Ne', 'Ar', 'S', 'Cl', 'Fe'}, 10);

    properties
        logAbundances
        elements
        elementRefence = 'H'
    end

    properties (Access = private)
        indexElementRefence
    end

    methods

        function obj = SolarAbundances(varargin)
            % Class constructor

            % Definitions
            defaultFilename = 'abundances.txt';
            
            % Create input parser
            ip = inputParser;
            addOptional(ip, 'filename', defaultFilename, @ischar);
            parse(ip, varargin{:});

            % Read abundances from filename
            [obj.logAbundances, obj.elements] = obj.read_abundances(ip.Results.filename);
            
            % Get index of the reference element
            obj.indexElementRefence = find_ind(obj.elements, obj.elementRefence);
        end

        function moles = abundances2moles(obj, elements, varargin)
            % Read solar abundances in log 10 scale and compute the initial molar
            % fractions in the mixture [-]
            %
            % Args:
            %     obj (SolarAbundances): SolarAbundances object
            %     elements (cell): List with the given elements
            %     filename (file): Filename with the data
            %
            % Optional Args:
            %     metallicity (float): Metallicity
            %
            % Returns:
            %     moles (float): moles relative to H of the remaining elements in the mixture
            %
            % Examples:
            %     * moles = SolarAbundances.abundances2moles({'H', 'He', 'C', 'N', 'O', 'Ne', 'Ar', 'S', 'Cl', 'Fe'})
            %     * moles = SolarAbundances.abundances2moles({'H', 'He', 'C', 'N', 'O', 'Ne', 'Ar', 'S', 'Cl', 'Fe'}, 10)
        
            % Unpack
            metallicity = unpack(varargin);

            % Get abundances from filename assuming unity metallicity
            logAbundances = obj.logAbundances;
            elementsDB = obj.elements;

            % Recompute with metallicity. NOTE: H and He do not change their abundances
            if nargin > 2
                index_H = find_ind(elementsDB, 'H');
                index_He = find_ind(elementsDB, 'He');
                index_change = 1:1:length(elementsDB);
                index_change([index_H, index_He]) = [];
                logAbundances(index_change) = logAbundances(index_change) + log10(metallicity);
            end
        
            % Reorganize abundances as the given element cell
            for i = length(elements):-1:1
                index(i) = find_ind(elementsDB, elements{i});
            end
        
            % Compute moles relative to H of the remaining elements in the mixture
            moles = 10.^(logAbundances(index) - logAbundances(obj.indexElementRefence));

            % SUB-PASS FUNCTIONS
            function metallicity = unpack(variable)
                % Unpack extra inputs
                metallicity = 1;
            
                if ~isempty(variable)
                    metallicity = variable{1};
                end
            
            end

        end

    end

    methods (Access = public, Static)
        
        function [abundances, elements] = read_abundances(filename)
            % Read solar abundances file
            %
            % Format: [number element, element, abundance, name, molar mass (g/mol)]
            %
            % Args:
            %     filename (file): Filename with the data
            %
            % Returns:
            %     Tuple containing
            %
            %     * abundances (float): Vector with the logarithmic base 10 solar abundances
            %     * elements (cell): List with the given elements
        
            % Open file
            fileID = fopen(filename, 'r');

            % Define the format of the data to read and the shape of the output array
            formatSpec = '%d %s %f %s %f';

            % Read the file data, filling output array, A, in column order. fscanf
            % reuses the format, formatSpec, throughout the file.
            A = textscan(fileID, formatSpec);

            % Close file
            fclose(fileID);

            % Get inputs
            abundances = A{1, 3};
            elements = A{1, 2};
        end

    end

end
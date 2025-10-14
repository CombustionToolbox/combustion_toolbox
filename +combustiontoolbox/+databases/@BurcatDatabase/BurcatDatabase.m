classdef BurcatDatabase < combustiontoolbox.databases.Database & handle
    % The :mat:func:`BurcatDatabase` class is used to store thermodynamic data from Burcat's database
    % using NASA's 9 coefficient polynomial fits.
    %
    % The :mat:func:`BurcatDatabase` object can be initialized as follows:
    %
    %       database = BurcatDatabase()
    %
    % This creates an instance of the :mat:func:`BurcatDatabase` class and initializes it with the
    % chemical species contained in Burcat's database.
    %
    % See also: :mat:func:`NasaDatabase`, :mat:func:`Database`

    methods (Access = public)
        
        function obj = BurcatDatabase(varargin)
            % Constructor 

            % Call superclass constructor
            obj@combustiontoolbox.databases.Database('name', 'Burcat', 'temperatureReference', 298.15, varargin{:});
        end

    end

    methods (Access = public, Static)

        function thermoMillennium2thermoNASA9(filenameInput, varargin)
            % Read Extended Third Millennium Thermodynamic Database of New NASA
            % Polynomials with Active Thermochemical Tables update and write a new
            % file compatible with thermo NASA 9 format
            %
            % Args:
            %     filenameInput (char): Filename of the thermoMillennium data
            %     filenameOutput (char): Filename of the thermoMillennium data as NASA9
            %
            % Optional args:
            %     * outDir (char): Output directory (default: 'databases')
            %     * suffix (char): Suffix to add to species names (default: '_M')
            %
            % Examples:
            %     * BurcatDatabase.thermoMillennium2thermoNASA9('thermo_millennium.inp')
            %     * BurcatDatabase.thermoMillennium2thermoNASA9('thermo_millennium.inp', 'filenameOutput', 'thermo_millennium_2_thermoNASA9.inp', 'outDir', 'databases', 'suffix', '_M')
            
            % CONSTANTS
            % MAX_CHAR = 80;

            % Default values
            defaultFilenameOutput = 'thermo_millennium_2_thermoNASA9.inp';
            defaultOutDir = 'databases';
            defaultSuffix = '_M';
            
            % Input parser
            p = inputParser;
            addRequired(p, 'filenameInput', @(x) ischar(x) || isstring(x));
            addOptional(p, 'filenameOutput', defaultFilenameOutput, @(x) ischar(x) || isstring(x));
            addParameter(p, 'outDir', defaultOutDir, @(x) ischar(x) || isstring(x));
            addParameter(p, 'suffix', defaultSuffix, @(x) ischar(x) || isstring(x));
            parse(p, filenameInput, varargin{:});

            % Set variables
            filenameInput = p.Results.filenameInput;
            filenameOutput = p.Results.filenameOutput;
            outDir = p.Results.outDir;
            SUFFIX = p.Results.suffix;

            % Definitions
            outPath = fullfile(outDir, filenameOutput);

            % Ensure output directory exists
            if ~exist(outDir, 'dir') && ~isempty(outDir)
                mkdir(outDir);
            end

            % Open source file
            [fid, msg] = fopen(filenameInput, 'r');
            if fid == -1
                error('Could not open input file: %s\nReason: %s', filenameInput, msg);
            end

            % Open destination file
            [fid_new, msg] = fopen(outPath, 'w');
            if fid_new == -1
                fclose(fid);
                error('Could not open output file: %s\nReason: %s', outPath, msg);
            end

            % Initialization
            tline = 1;
            FLAG_NEW_SPECIES = true;
            N_SUFFIX = length(SUFFIX);
            speciesMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
        
            while tline ~= -1
                tline = fgetl(fid);
        
                if ~ischar(tline)
                    break
                elseif isempty(tline)
                    FLAG_NEW_SPECIES = true;
                    tline = 1;
                    continue
                elseif tline(1) == '!'
                    continue
                elseif contains(tline, 'see')
                    continue
                elseif ~isempty(regexp(tline, '^[A-Za-z][^\s]*\s+', 'once'))
                    FLAG_NEW_SPECIES = true;
                end

                % Extract species data
                if FLAG_NEW_SPECIES
                    ind_space = regexp(tline, '\s');
                    ind_next = regexp(tline(ind_space(1):end), '\S');
                    N = 18 - ind_space(1) - 1 - N_SUFFIX;
                    white_spaces = blanks(N);
                    species = tline(1:ind_space(1) - 1);

                    if contains(tline, ' cr ', 'IgnoreCase', true)
                        species = strcat(species, '(cr)');
                    elseif contains(tline, ' liq ', 'IgnoreCase', false) && ~contains(species, '(l)', 'IgnoreCase', true)
                        species = strcat(species, '(L)');
                    end

                    if contains(tline, 'excited', 'IgnoreCase', true)
                        species = strcat(species, '(exc)');
                    elseif contains(tline, 'singlet', 'IgnoreCase', true)
                        species = strcat(species, '(slet)');
                    elseif contains(tline, 'doublet', 'IgnoreCase', true)
                        species = strcat(species, '(dlet)');
                    elseif contains(tline, 'triplet', 'IgnoreCase', true)
                        species = strcat(species, '(tlet)');
                    elseif contains(tline, 'quartet', 'IgnoreCase', true)
                        species = strcat(species, '(qtet)');
                    end
        
                    species = replace(species, '*', ' ');
                    species = replace(species, '=', '_');
                    species = replace(species, 'Al', 'AL');
                    species = replace(species, 'Cl', 'CL');
                    species = replace(species, 'Tl', 'TL');
                    species = replace(species, 'Fl', 'FL');
                    species = replace(species, '(liq)', 'liq');
                    species = replace(species, 'liq', '(liq)');
                    species = replace(species, '(l)', '(L)');
                    
                    % Check container
                    if isKey(speciesMap, species)
                        % Repeated species in the container
                        speciesMap(species) = speciesMap(species) + 1;
                        species = sprintf('%s_num%d', species, speciesMap(species));
                    else
                        % Add species to the container
                        speciesMap(species) = 1;
                    end
                    
                    % Add suffix
                    species = strcat(species, SUFFIX);
                    
                    fprintf(fid_new, '%s%s%s\n', species, white_spaces, tline(ind_next(1) + ind_space(1) - 1:end));
                    FLAG_NEW_SPECIES = false;
                else
                    fprintf(fid_new, '%s\n', tline);
                end
        
            end
        
            fclose(fid);
            fclose(fid_new);
        end
        
    end

end
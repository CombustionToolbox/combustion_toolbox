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

        function thermoMillennium_2_thermoNASA9(filenameInput, varargin)
            % Read Extended Third Millennium Thermodynamic Database of New NASA
            % Polynomials with Active Thermochemical Tables update and write a new
            % file compatible with thermo NASA 9 format
            %
            % Args:
            %     filenameInput (char): Filename of the thermoMillennium data
            %     filenameOutput (char): Filename of the thermoMillennium data as NASA9
            
            % Default values
            filenameOutput = 'thermo_millennium_2_thermoNASA9.inp';

            if nargin > 1
                filenameOutput = varargin{1};
            end

            
            fid = fopen(filenameInput, 'r');
            fid_new = fopen(strcat('databases/', filenameOutput), 'w');
            tline = 1;
            FLAG_NEW_SPECIES = true;
            MAX_CHAR = 80;
            SUFFIX = '_M';
            N_SUFFIX = length(SUFFIX);
        
            while tline ~= -1
                tline = fgetl(fid);
        
                if ~ischar(tline)
                    break
                elseif isempty(tline)
                    FLAG_NEW_SPECIES = true;
                    tline = 1;
                    continue
                end
        
                if tline(1) == '!'
                    continue
                end
        
                if FLAG_NEW_SPECIES
                    ind_space = regexp(tline, '\s');
                    ind_next = regexp(tline(ind_space(1):end), '\S');
                    N = 18 - ind_space(1) - 1 - N_SUFFIX;
                    white_spaces = blanks(N);
                    species = strcat(tline(1:ind_space(1) - 1), SUFFIX);
        
                    if contains(tline, 'excited', 'IgnoreCase')
                        species = replace(species, SUFFIX, strcat('bexb', SUFFIX));
                    end
        
                    if contains(tline, 'singlet', 'IgnoreCase')
                        species = replace(species, SUFFIX, strcat('bsingletb', SUFFIX));
                    elseif contains(tline, 'doublet', 'IgnoreCase')
                        species = replace(species, SUFFIX, strcat('bdoubletb', SUFFIX));
                    elseif contains(tline, 'triplet', 'IgnoreCase')
                        species = replace(species, SUFFIX, strcat('btripletb', SUFFIX));
                    elseif contains(tline, 'quartet', 'IgnoreCase')
                        species = replace(species, SUFFIX, strcat('bquartetb', SUFFIX));
                    end
        
                    species = replace(species, '*', ' ');
                    species = replace(species, '=', '_');
                    species = replace(species, 'Al', 'AL');
                    species = replace(species, 'Cl', 'CL');
                    species = replace(species, 'Tl', 'TL');
                    species = replace(species, 'Fl', 'FL');
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
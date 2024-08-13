classdef NasaDatabase < combustiontoolbox.databases.Database & handle
    % The :mat:func:`NasaDatabase` class is used to store thermodynamic data from NASA's database
    % using NASA's 9 coefficient polynomial fits.
    %
    % The :mat:func:`NasaDatabase` object can be initialized as follows:
    %
    %       database = NasaDatabase()
    %
    % This creates an instance of the :mat:func:`NasaDatabase` class and initializes it with the
    % chemical species contained in NASA's database.
    %
    % See also: :mat:func:`BurcatDatabase`, :mat:func:`Database`

    methods (Access = public)
        
        g0 = species_g0_NASA(obj, species, temperature)

        function obj = NasaDatabase(varargin)
            % Constructor
            %
            % Optional Args:
            %   varargin (optional): key-value pairs to initialize the database
            %
            % Returns:
            %     obj (NasaDatabase): Object with NASA's database
            %
            % Examples:
            %     * db = combustiontoolbox.databases.NasaDatabase();
            %     * db = combustiontoolbox.databases.NasaDatabase('filename', 'DB.mat');

            % Call superclass constructor
            obj@combustiontoolbox.databases.Database('name', 'NASA', 'temperatureReference', 298.15, varargin{:});
        end

        function DB_master = generateDatabaseMaster(obj)
            % Generate Master Database (DB_master) with the thermodynamic
            % data of the chemical species
            %
            % Args:
            %     obj (NasaDatabase): NasaDatabase object
            %
            % Returns:
            %     DB_master (struct): Database with the thermodynamic data of the chemical species
            %
            % Example:
            %     * DB_master = generateDatabaseMaster('thermo_CT.inp')
            
            % Generate master database
            DB_master = obj.getDatabaseMaster(obj.thermoFile);
        
            fprintf('OK!\n');
        end

        function obj = generateDatabase(obj, varargin)
            % Generate Database with thermochemical interpolation curves
            % from the data extracted from the thermoFile
            % 
            % Args:
            %     obj (NasaDatabase): NasaDatabase object
            %
            % Optional Args:
            %     * listSpecies (cell): List of species to generate the database
            %
            % Returns:
            %     obj (NasaDatabase): NasaDatabase object with thermochemical interpolation curves
            %
            % Examples:
            %     * DB = generateDatabase(NasaDatabase());
            %     * DB = generateDatabase(NasaDatabase(), {'N2', 'O2', 'NO', 'O', 'N'});
            
            % Get master database from the thermoFile
            DB_master = getDatabaseMaster(obj, obj.thermoFile);

            % Unpack inputs
            if nargin > 1
                listSpecies = varargin{1};

                % Remove not required species
                DB_master = rmfield(DB_master, setdiff(fieldnames(DB_master), listSpecies));
            else
                listSpecies = fieldnames(DB_master);
            end

            % Definitions
            numSpecies = length(listSpecies);

            % Control message
            fprintf('Generating %s database with thermo ... ', obj.name);
            
            % Compute interpolation curves for each species
            for i = 1:numSpecies
                species = obj.fullname2name(listSpecies{i});

                if ~isfield(DB_master, species)
                    fprintf(['\n- Species ''', listSpecies{i}, ''' does not exist as a field in species structure ... ']);
                    continue
                end
                
                % Initialization
                temp = DB_master.(species);

                % Get data
                Tintervals = DB_master.(species).Tintervals;
                Trange = DB_master.(species).Trange;
                
                % Get thermodynamic data from the species that cannot be evaluated at different temperatures
                if temp.Tintervals == 0
                    Tref = Trange(1);

                    [Cp0, Hf0, H0, Ef0, S0, DfG0] = obj.getSpeciesThermo(obj, DB_master, listSpecies{i}, Tref, obj.units);
                    
                    temp.hf = Hf0;
                    temp.ef = Ef0;
                    temp.Tref = Tref;
                    temp.T = Tref;

                    % Interpolation curves (constant values)
                    temp.cpcurve = griddedInterpolant([Tref, Tref + 1], [Cp0, Cp0], 'linear', 'linear');
                    temp.h0curve = griddedInterpolant([Tref, Tref + 1], [H0, H0], 'linear', 'linear');
                    temp.s0curve = griddedInterpolant([Tref, Tref + 1], [S0, S0], 'linear', 'linear');
                    temp.g0curve = griddedInterpolant([Tref, Tref + 1], [DfG0, DfG0], 'linear', 'linear');

                    % Store the species data in the SpeciesDB
                    DB.(species) = temp;

                    continue
                end
                
                % Get thermodynamic data from the species that can be evaluated at different temperatures
                [~, Hf0, ~, Ef0, ~, ~] = obj.getSpeciesThermo(obj, DB_master, listSpecies{i}, obj.temperatureReference, obj.units);
                
                temp.hf = Hf0;
                temp.ef = Ef0;
                temp.Tref = obj.temperatureReference;

                Tmin = Trange{1}(1);
                Tmax = Trange{Tintervals}(2);
                T_vector = linspace(Tmin, Tmax, obj.pointsTemperature);

                [Cp0_vector, ~, H0_vector, ~, S0_vector, ~] = obj.getSpeciesThermo(obj, DB_master, listSpecies{i}, T_vector, obj.units);
                DfG0_vector = H0_vector - T_vector .* S0_vector;

                temp.T = T_vector;

                % Interpolation curves
                temp.cpcurve = griddedInterpolant(T_vector, Cp0_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.h0curve = griddedInterpolant(T_vector, H0_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.s0curve = griddedInterpolant(T_vector, S0_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.g0curve = griddedInterpolant(T_vector, DfG0_vector, obj.interpolationMethod, obj.extrapolationMethod);

                % Coefficients NASA's 9 polynomial fits
                temp.Tintervals = DB_master.(species).Tintervals;
                temp.Trange = DB_master.(species).Trange;
                temp.Texponents = DB_master.(species).Texponents;
                temp.a = DB_master.(species).a;
                temp.b = DB_master.(species).b;

                % Store the species data in the database
                DB.(species) = temp;
            end
            
            % Assign data
            obj.species = DB;
        end
        
        function [cp, cv, h0, DhT, e0, DeT, s0, g0] = getSpeciesThermoFull(obj, DB, species, temperature, units)
            % Compute thermodynamic function using NASA's 9 polynomials
            %
            % Args:
            %     species (char): Chemical species
            %     temperature (float): Range of temperatures to evaluate [K]
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %     units (char): Label indicating mass [kg] or molar [mol] units
            %
            % Returns:
            %     Tuple containing
            %
            %     * cp  (float): Specific heat at constant pressure in units basis [J/(units-K)]
            %     * cv  (float): Specific heat at constant volume in units basis   [J/(units-K)]
            %     * h0  (float): Enthalpy in units basis [J/units]
            %     * DhT (float): Thermal enthalpy in units basis [J/units]
            %     * e0  (float): Internal energy in units basis [J/units]
            %     * DeT (float): Thermal internal energy in units basis [J/units]
            %     * s0  (float): Entropy in units basis [J/(units-K)]
            %     * g0  (float): Gibbs energy in units basis [J/units]
            %
            % Example:
            %     [cp, cv, h0, DhT, e0, DeT, s0, g0] = speciesThermoFull('H2O', 300:100:6000, DB)
            
            % Definitions
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]

            % Compute thermodynamic functions included in the propertiesMatrix
            [cp, hf, h0, ef, s0, g0] = obj.getSpeciesThermo(obj, DB, species, temperature, units);
            
            % Unpack NASA's polynomials coefficients
            [~, ~, Trange, ~, Tintervals, ~, ~, ~, ~] = obj.getCoefficients(species, DB);

            % Compute additional thermodynamic function not included in getSpeciesThermo
            if Tintervals > 0
                e0 = (ef + (h0 - hf) - (1 - phase) * R0 * (T - obj.temperatureReference));
                cv = cp - R0;
                DhT = h0 - hf;
                DeT = e0 - ef;
            else
                Tref = Trange(1);
                cv = zeros(1, N);
                e0 = hf - Delta_n * R0 * Tref;
                s0 = zeros(1, N);
                DhT = zeros(1, N);
                DeT = zeros(1, N);
            end
        
            % Change units
            % h0 = h0 * 1e-3;   % [kJ/mol]
            % DhT = DhT * 1e-3; % [kJ/mol]
            % e0 = e0 * 1e-3;   % [kJ/mol]
            % DeT = DeT * 1e-3; % [kJ/mol]
            % g0 = g0 * 1e-3;   % [kJ/mol]
            % s0 = s0 * 1e-3;   % [kJ/(mol-K)]
        end

    end

    methods (Access = public, Static)
        
        function name = fullname2name(species)
            % Get full name of the given species
            %
            % Args:
            %     species (char): Chemical species
            %
            % Returns:
            %     name (char): Full name of the given species
        
            FLAG_MILLENIUM = false;
        
            if contains(species, '_M')
                species = strrep(species, '_M', '');
                FLAG_MILLENIUM = true;
            end
        
            name = species;
        
            if isempty(name)
                return
            end
        
            if name(end) == '+'
                name = [name(1:end - 1) 'plus'];
            elseif name(end) == '-'
                name = [name(1:end - 1) 'minus'];
            end
        
            ind = regexp(name, '[()]');
            name(ind) = 'b';
            ind = regexp(name, '[.,+-]');
            name(ind) = '_';
        
            if regexp(name(1), '[0-9]')
                name = ['num_' name];
            end
        
            ind = regexp(name, '\x27');
            name(ind) = '_';
        
            if FLAG_MILLENIUM
                name = strcat(name, '_M');
            end
        
        end

        function [a, b, Trange, Texponents, Tintervals, phase, hf0, W, FLAG_REFERENCE] = getCoefficients(species, DB)
            % Unpack NASA's polynomials coefficients from database
            %
            % Args:
            %     species (char): Chemical species
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %
            % Returns:
            %     Tuple containing
            %
            %     * a (cell): Temperature coefficients
            %     * b (cell): Integration constants
            %     * Trange (cell): Ranges of temperatures [K]
            %     * Texponents (cell): Exponent coefficients
            %     * Tintervals (float): Number of intervals of temperatures
            %     * phase (float): 0 or 1 indicating gas or condensed phase, respectively
            %     * hf0 (float): Enthalpy of formation [J/mol]
            %     * W (float): Molecular weight [kg/mol]
            %     * FLAG_REFERENCE (bool): Flag indicating species is a reference element/species
            %
            % Example:
            %     [a, b, Trange, Texponents, Tintervals, phase, hf0, W, FLAG_REFERENCE] = getCoefficients('H2O', DB)
        
            a = DB.(species).a;
            b = DB.(species).b;
            Trange = DB.(species).Trange;
            Texponents = DB.(species).Texponents;
            Tintervals = DB.(species).Tintervals;
            phase = DB.(species).phase;
            hf0 = DB.(species).hf;
            W = DB.(species).W;
            FLAG_REFERENCE = DB.(species).FLAG_REFERENCE;
        end

    end
    
    methods (Access = private)

        function DB_master = getDatabaseMaster(obj, thermoFile)
            % Generate Master Database (DB_master) with the thermodynamic 
            % data of the chemical species
            %
            % Args:
            %     thermoFile (char): path of thermoFile
            %
            % Returns:
            %     DB_master (struct): Database with the thermodynamic data of the chemical species
    
            % Load database
            fid = fopen(thermoFile); 
    
            switch thermoFile
                case 'thermo_CT.inp'
                    msg = 'Loading NASA database ... ';
                otherwise
                    msg = 'Loading an unkown database ... ';
            end
        
            fprintf(msg)
            line = 0;
        
            while line < 1e4
                tline = fgetl(fid);
        
                if ~ischar(tline)
                    break
                end
        
                if isempty(regexp(tline, '\S', 'once'))
                    continue
                end
        
                if tline(1) == '!'
                    continue
                end
        
                if contains(tline, 'thermo')
                    tline = fgetl(fid);
                    continue
                end
        
                if contains(tline, 'END')
                    continue
                end
        
                line = line + 1;
                temp = combustiontoolbox.core.Species();
                temp.fullname = sscanf(tline(1:16), '%s');
                temp.name = obj.fullname2name(temp.fullname);
                temp.comments = tline(19:end);

                if contains(temp.comments, 'Ref-')
                    temp.FLAG_REFERENCE = true;
                end

                tline = fgetl(fid);
                temp.Tintervals = str2double(tline(1:2));
                temp.refCode = tline(4:9);
                temp.formula = tline(11:50);
                temp.phase = str2double(tline(51:52));
                temp.W = str2double(tline(53:65)) * 1e-3; % [kg/mol]
                temp.hf = str2double(tline(66:80));
        
                if temp.Tintervals == 0
                    tline = fgetl(fid);
                    temp.Trange = str2num(tline(1:22)); %#ok<ST2NM>
                    temp.Texponents = str2num(tline(24:63)); %#ok<ST2NM>
                    temp.hftoh0 = str2double(tline(66:end));
                end
        
                for Tinterval = 1:temp.Tintervals
                    tline = fgetl(fid);
                    temp.Trange{Tinterval} = str2num(tline(1:22)); %#ok<ST2NM>
                    temp.Texponents{Tinterval} = str2num(tline(24:63)); %#ok<ST2NM>
                    temp.hftoh0{Tinterval} = str2double(tline(66:end));
        
                    tline = fgetl(fid);
                    a1 = str2num(tline(1:16)); %#ok<ST2NM>
                    a2 = str2num(tline((1:16) + 16)); %#ok<ST2NM>
                    a3 = str2num(tline((1:16) + 32)); %#ok<ST2NM>
                    a4 = str2num(tline((1:16) + 48)); %#ok<ST2NM>
                    a5 = str2num(tline((1:16) + 64)); %#ok<ST2NM>
        
                    tline = fgetl(fid);
                    a6 = str2num(tline(1:16)); %#ok<ST2NM>
                    a7 = str2num(tline((1:16) + 16)); %#ok<ST2NM>
                    a8 = 0; %str2num(tline((1:16)+32));
                    b1 = str2num(tline((1:16) + 48)); %#ok<ST2NM>
                    b2 = str2num(tline((1:16) + 64)); %#ok<ST2NM>
                    temp.a{Tinterval} = [a1 a2 a3 a4 a5 a6 a7 a8];
                    temp.b{Tinterval} = [b1 b2];
                end
        
                DB_master.(temp.name) = temp;
            end
        
            fclose(fid);
        end

    end

    methods (Access = private, Static)

        function [cp0, hf, h0, ef, s0, g0] = getSpeciesThermo(obj, DB, species, temperature, units)
            % Calculates the thermodynamic properties of any species included in the NASA database
            %
            % Args:
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %     species (char): Chemical species
            %     T (float): Temperature [K]
            %     units (char): Label indicating mass [kg] or molar [mol] units
            %
            % Returns:
            %     Tuple containing
            %
            %     * cP0 (float): Specific heat at constant pressure [J/(mol-k)]
            %     * hf0 (float): Enthalpy of formation [J/mol]
            %     * h0 (float):  Enthalpy [J/mol]
            %     * ef0 (float): Internal energy of formation [J/mol]
            %     * s0 (float):  Entropy [J/(mol-k)]
            %     * Dg0 (float): Gibbs energy [J/mol]
            %
            % Example:
            %     * [formula, mm, cP0, hf0, h0, ef0, s0, g0] = getSpeciesThermo(obj, DB, 'CO', 1000, 'molar')
            %
            %     * formula = 'C   1.00O   1.00    0.00    0.00    0.00'
            %     * mm  = 28.0101
            %     * cP0 = 33.1788
            %     * hf0 = -1.1054e+05
            %     * h0  = -8.8848e+04
            %     * ef0 = -1.1177e+05
            %     * s0  = 234.5409
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Definitions
            R0 = combustiontoolbox.common.Constants.R0; % Universal gas constant [J/(K mol)]
            N = length(temperature);

            % Unpack NASA's polynomials coefficients
            [a, b, Trange, Texponents, Tintervals, phase, hf, W, FLAG_REFERENCE] = obj.getCoefficients(species, DB);
        
            % Get elements
            elements = combustiontoolbox.core.Elements();

            % Get element matrix of the species
            elementMatrix = DB.(species).getElementMatrix(elements.listElements);
        
            % Compute change in moles of gases during the formation reaction of a
            % mole of that species starting from the elements in their reference state
            Delta_n = obj.getChangeMolesGasReaction(elements, elementMatrix, phase);
        
            % If temperature interval is zero the species is only a reactant. In that case,
            % determine it's reference temperature Tref. For noncryogenic reactants, assigned
            % enthalpies are given at 298.15 K. For cryogenic liquids, assigned enthalpies
            % are given at their boiling points instead of 298.15 K
            if Tintervals == 0
                Tref = Trange(1);
                cp0 = zeros(1, N);
                h0 = hf * ones(1, N);
                ef = hf - Delta_n * R0 * Tref;
                s0 = zeros(1, N);
                g0 = h0;
        
                if strcmpi(units, 'mass')
                    hf = molar2mass(hf, W);
                    ef = molar2mass(ef, W);
                end

                return
            end

            % Compute thermodynamic properties
            for i = N:-1:1
                T = temperature(i);

                if Tintervals > 0
                    Tref = obj.temperatureReference;
            
                    % Get temperature interval
                    Tinterval = obj.getIndexTempereratureInterval(species, T, DB);

                    % Compute thermodynamic function
                    cp0(i) = R0 * sum(a{Tinterval} .* T.^Texponents{Tinterval});
                    h0(i) = R0 * T * (sum(a{Tinterval} .* T.^Texponents{Tinterval} .* [-1 log(T) 1 1/2 1/3 1/4 1/5 0]) + b{Tinterval}(1) / T);
                    s0(i) = R0 * (sum(a{Tinterval} .* T.^Texponents{Tinterval} .* [-1/2 -1 log(T) 1 1/2 1/3 1/4 0]) + b{Tinterval}(2));
                    ef(i) = hf - Delta_n * R0 * Tref;
            
                    % Compute the standar gibbs free energy of formation at the specified
                    % temperature. This enforces us to consider explicitely the formation
                    % reaction from the elements in their reference states at room
                    % temperature, unless the species is precisely an element in its
                    % reference state, in which case the standard gibbs free energy of
                    % formation is identically zero.
                    if ~FLAG_REFERENCE
                        g0(i) = h0(i) - T .* s0(i);
                    else
                        g0(i) = 0;
                    end
            
                end

            end
            
            if strcmpi(units, 'mass')
                cp0 = molar2mass(cp0, W);
                hf = molar2mass(hf, W);
                ef = molar2mass(ef, W);
                h0 = molar2mass(h0, W);
                s0 = molar2mass(s0, W);
                g0 = molar2mass(g0, W);
            end

            % SUB-PASS FUNCTIONS
            function value = molar2mass(value, W)
                % Change molar units [mol] to mass units [kg]
                value = value / W;
            end

        end

        function Delta_n = getChangeMolesGasReaction(elements, elementMatrix, phase)
            % In order to compute the internal energy of formation from the enthalpy of
            % formation of a given species, we must determine the change in moles of
            % gases during the formation reaction of a mole of that species starting
            % from the elements in their reference state. 
            % 
            % Notes:
            %     The only elements that are stable as diatomic gases are elements
            %     1 (H), 8 (N), 9 (O), 10 (F), and 18 (Cl). The remaining elements that
            %     are stable as (monoatomic) gases are the noble gases He (3), Ne (11),
            %     Ar (19), Kr (37), Xe (55), and Rn (87), which do not form any compound.
            %
            % Args:
            %     elements (Elements): Elements object
            %     elementMatrix (float): Element matrix of the species
            %     phase (float): 0 or 1 indicating gas or condensed species
            %
            % Optional Args:
            %     * elements (Elements): Elements object
            %
            % Returns:
            %     Delta_n (float): Change in moles of gases during the formation reaction of a mole of that species starting from the elements in their reference state
            %
            % Example:
            %     Delta_n = getChangeMolesGasReaction(elements, element_matrix, phase)

            Delta_n_per_mole = sum(elementMatrix(1,:) == [elements.indexH, elements.indexN, elements.indexO, elements.indexF, elements.indexCl]') / 2 ... 
                             + sum(elementMatrix(1,:) == [elements.indexHe, elements.indexNe, elements.indexAr, elements.indexKr, elements.indexXe, elements.indexRn]');
            Delta_n = 1 - phase - dot(Delta_n_per_mole, elementMatrix(2,:));
        end

        function Tinterval = getIndexTempereratureInterval(species, T, DB)
            % Get interval of the NASA's polynomials from the Database (DB) for the
            % given species and temperature [K].
            %
            % Args:
            %     species (char): Chemical species
            %     T (float): Temperature [K]
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %
            % Returns:
            %     Tinterval (float): Index of the interval of temperatures
            
            for i = 1:DB.(species).Tintervals
        
                if (T >= DB.(species).Trange{i}(1)) && (T <= DB.(species).Trange{i}(2))
                    break
                end
        
            end
        
            Tinterval = i;
        end

    end

end
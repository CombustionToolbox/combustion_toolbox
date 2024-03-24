classdef Database

    properties
        name
        species
        filename
        filenameMaster
        interpolationMethod
        extrapolationMethod 
        units
        pointsTemperature
        temperatureReference
        thermoFile
    end

    properties (Hidden)
        id
    end

    properties (Dependent)
        listSpecies
        numSpecies
    end
    
    methods

        function value = get.listSpecies(obj)
            value = fieldnames(obj.species);
        end

        function value = get.numSpecies(obj)
            value = length(obj.listSpecies);
        end

    end

    methods (Access = public)
        
        function obj = Database(varargin)
            % Constructor
            
            % Definitions
            defaultName = 'Database';
            defaultFilename = 'DB.mat';
            defaultFilenameMaster = 'DB_master.mat';
            defaultInterpolationMethod = 'pchip';
            defaultExtrapolationMethod = 'linear';
            defaultUnits = 'molar';
            defaultPointsTemperature = 200;
            defaultTemperatureReference = 298.15; % [K]
            defaultThermoFile = 'thermo_CT.inp';

            % Parse function inputs
            ip = inputParser;
            addParameter(ip, 'name', defaultName, @ischar);
            addParameter(ip, 'species', [], @iscell);
            addParameter(ip, 'filename', defaultFilename, @ischar);
            addParameter(ip, 'filenameMaster', defaultFilenameMaster, @ischar);
            addParameter(ip, 'interpolationMethod', defaultInterpolationMethod, @ischar);
            addParameter(ip, 'extrapolationMethod', defaultExtrapolationMethod, @ischar);
            addParameter(ip, 'units', defaultUnits, @ischar);
            addParameter(ip, 'pointsTemperature', defaultPointsTemperature, @isnumeric);
            addParameter(ip, 'temperatureReference', defaultTemperatureReference, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
            addParameter(ip, 'thermoFile', defaultThermoFile);
            parse(ip, varargin{:});
            
            % Set properties
            obj.name = ip.Results.name;
            obj.species = ip.Results.species;
            obj.filename = ip.Results.filename;
            obj.filenameMaster = ip.Results.filenameMaster;
            obj.interpolationMethod = ip.Results.interpolationMethod;
            obj.extrapolationMethod = ip.Results.extrapolationMethod;
            obj.units = ip.Results.units;
            obj.pointsTemperature = ip.Results.pointsTemperature;
            obj.temperatureReference = ip.Results.temperatureReference;
            obj.thermoFile = ip.Results.thermoFile;

            % Generate ID
            obj = obj.generate_id();

            % Check if database is in cached and the id matches
            persistent cachedDatabase;
            if ~isempty(cachedDatabase) && isequal(cachedDatabase.id, obj.id)
                obj = cachedDatabase;
            else
                % Load database
                obj = obj.load();
                % Cache database
                cachedDatabase = obj;
            end

        end

        function obj = load(obj, varargin)
            %

            if nargin > 1
                obj.filename = varargin{1};
            end

            % Load database
            if exist(obj.filename, 'file')
                fprintf('%s database with thermo loaded from the main path ... ', obj.name);
                load(obj.filename, 'DB');
                obj = DB;
            else
                DB_master = obj.generateDatabaseMaster();
                DB = obj.generate_DB(DB_master);
                obj.species = DB;
            end

            % Status
            fprintf('OK!\n');
        end

        function save(obj)
            % Save database
            fprintf('Work in progress\n');
        end
        
        function DB_master = generateDatabaseMaster(obj)
            % Generate Master Database (DB_master) with the thermodynamic data of
            % the chemical species
            %
            % Args:
            %     obj (Database): 
            %
            % Returns:
            %     DB_master (struct): Database with the thermodynamic data of the chemical species
            %
            % Example:
            %     * DB_master = generate_DB_master('thermo_CT.inp')
            
            % Load master database
            if exist(obj.filenameMaster, 'file')
                fprintf('Loading NASA database ... ')
                load('DB_master.mat', 'DB_master');
            else
                DB_master = getDatabaseMaster(obj.thermoFile);        
            end
        
            fprintf('OK!\n');

            % SUB-PASS FUNCTIONS
            function DB_master = getDatabaseMaster(thermoFile)
                % 

                % Import packages
                import combustiontoolbox.core.Species
                
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
                    temp = Species();
                    temp.fullname = sscanf(tline(1:16), '%s');
                    temp.name = obj.FullName2name(temp.fullname);
                    temp.comments = tline(19:end);
                    tline = fgetl(fid);
                    temp.Tintervals = str2double(tline(1:2));
                    temp.refCode = tline(4:9);
                    temp.formula = tline(11:50);
                    temp.phase = str2double(tline(51:52));
                    temp.W = str2double(tline(53:65));
                    temp.hf = str2double(tline(66:80));
            
                    if temp.Tintervals == 0
                        tline = fgetl(fid);
                        temp.Trange = str2num(tline(1:22)); %#ok<ST2NM>
                        temp.Texponents = str2num(tline(24:63)); %#ok<ST2NM>
                        temp.hftoh0 = str2double(tline(66:end));
                    end
            
                    for ctInterval = 1:temp.Tintervals
                        tline = fgetl(fid);
                        temp.Trange{ctInterval} = str2num(tline(1:22)); %#ok<ST2NM>
                        temp.Texponents{ctInterval} = str2num(tline(24:63)); %#ok<ST2NM>
                        temp.hftoh0{ctInterval} = str2double(tline(66:end));
            
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
                        temp.a{ctInterval} = [a1 a2 a3 a4 a5 a6 a7 a8];
                        temp.b{ctInterval} = [b1 b2];
                    end
            
                    DB_master.(temp.name) = temp;
                end
            
                fclose(fid);
            end
            
        end

        function DB = generate_DB(obj, DB_master)
            % Generate Database with thermochemical interpolation curves for
            % the species contained in DB_master
            
            % Definitions
            LS = fieldnames(DB_master);
            NS = length(LS);

            % Control message
            fprintf('Generating %s database with thermo ... ', obj.name);
            
            % Compute interpolation curves for each species
            for i = 1:NS
                species = obj.FullName2name(LS{i});

                if ~isfield(DB_master, species)
                    fprintf(['\n- Species ''', LS{i}, ''' does not exist as a field in DB_master structure ... ']);
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

                    [Cp0, Hf0, H0, Ef0, S0, DfG0] = obj.getSpeciesProperties(obj, DB_master, LS{i}, Tref, obj.units, 0);
                    
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
                [~, Hf0, ~, Ef0, ~, ~] = obj.getSpeciesProperties(obj, DB_master, LS{i}, obj.temperatureReference, obj.units, 0);
                
                temp.hf = Hf0;
                temp.ef = Ef0;
                temp.Tref = obj.temperatureReference;

                Tmin = Trange{1}(1);
                Tmax = Trange{Tintervals}(2);
                T_vector = linspace(Tmin, Tmax, obj.pointsTemperature);

                for j = obj.pointsTemperature:-1:1
                    [Cp0, ~, H0, ~, S0, ~] = obj.getSpeciesProperties(obj, DB_master, LS{i}, T_vector(j), obj.units, 0);
                    cp_vector(j) = Cp0;
                    h0_vector(j) = H0;
                    s0_vector(j) = S0;
                    g0_vector(j) = H0 - T_vector(j) * S0;
                end

                temp.T = T_vector;

                % Interpolation curves
                temp.cpcurve = griddedInterpolant(T_vector, cp_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.h0curve = griddedInterpolant(T_vector, h0_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.s0curve = griddedInterpolant(T_vector, s0_vector, obj.interpolationMethod, obj.extrapolationMethod);
                temp.g0curve = griddedInterpolant(T_vector, g0_vector, obj.interpolationMethod, obj.extrapolationMethod);

                % Coefficients NASA's 9 polynomial fits
                temp.Tintervals = DB_master.(species).Tintervals;
                temp.Trange = DB_master.(species).Trange;
                temp.Texponents = DB_master.(species).Texponents;
                temp.a = DB_master.(species).a;
                temp.b = DB_master.(species).b;

                % Store the species data in the database
                DB.(species) = temp;
                
            end

        end
    
    end
    
    methods (Access = protected)

        function obj = generate_id(obj)
            % Concatenate input arguments to create a unique identifier string
            value =  [obj.name, num2str(obj.species), obj.filename, obj.filenameMaster, ...
                      obj.interpolationMethod, obj.extrapolationMethod, obj.units, ...
                      num2str(obj.pointsTemperature), num2str(obj.temperatureReference),...
                      obj.thermoFile];

            obj.id = combustiontoolbox.utils.generate_id(value);
        end

    end

    methods (Access = private, Static)
        
        function name = FullName2name(species)
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

        function elementsTemperatureReference = setElementsTemperatureReference()
        
            elementsTemperatureReference = {
                        'Ar [200-20000]';
                        'CL2 [200-6000]';
                        'F2 [200-6000]';
                        'H2 [200-20000]';
                        'He [200-20000]';
                        'Kr [200-20000]';
                        'N2 [200-20000]';
                        'Ne [200-20000]';
                        'O2 [200-20000]';
                        'Rn [200-20000]';
                        'Xe [200-20000]';
                        'Ag(cr) [200-1235.08]';
                        'Ag(L) [1235.08-6000]';
                        'AL(cr) [200-933.61]';
                        'AL(L) [933.61-6000]';
                        'B(b) [200-2350]';
                        'B(L) [2350-6000]';
                        'Ba(cr) [80-1000]';
                        'Ba(L) [1000-6000]';
                        'Be(a) [100-1543]';
                        'Be(b) [1543-1563]';
                        'Be(L) [1563-6000]';
                        'Br2(L) [265.9-6000]';
                        'C(gr) [200-6000]';
                        'Ca(a) [200-716]';
                        'Ca(b) [716-1115]';
                        'Ca(L) [1115-6000]';
                        'Cd(cr) [100-594.258]';
                        'Cd(L) [594.258-6000]';
                        'Co(a) [200-700.1]';
                        'Co(b) [700.1-1394]';
                        'Co(b) [1394-1768]';
                        'Co(L) [1768-6000]';
                        'Cr(cr) [200-311.5]';
                        'Cr(cr) [311.5-2130]';
                        'Cr(L) [2130-6000]';
                        'Cs(cr) [100-301.59]';
                        'Cs(L) [301.59-2000]';
                        'Cu(cr) [200-1358]';
                        'Cu(L) [1358-6000]';
                        'Fe(a) [200-1042]';
                        'Fe(a) [1042-1184]';
                        'Fe(c) [1184-1665]';
                        'Fe(d) [1665-1809]';
                        'Fe(L) [1809-6000]';
                        'Ga(cr) [100-302.92]';
                        'Ga(L) [302.92-6000]';
                        'Ge(cr) [200-1211.4]';
                        'Ge(L) [1211.4-6000]';
                        'Hg(cr) [100-234.29]';
                        'Hg(L) [234.29-2000]';
                        'I2(cr) [200-386.75]';
                        'I2(L) [386.75-6000]';
                        'In(cr) [100-429.784]';
                        'In(L) [429.784-6000]';
                        'K(cr) [200-336.86]';
                        'K(L) [336.86-2200]';
                        'Li(cr) [200-453.69]';
                        'Li(L) [453.69-6000]';
                        'Mg(cr) [100-923]';
                        'Mg(L) [923-6000]';
                        'Mn(a) [200-980]';
                        'Mn(b) [980-1361]';
                        'Mn(c) [1361-1412]';
                        'Mn(d) [1412-1519]';
                        'Mn(L) [1519-6000]';
                        'Mo(cr) [200-2896]';
                        'Mo(L) [2896-6000]';
                        'Na(cr) [200-371.01]';
                        'Na(L) [371.01-2300]';
                        'Nb(cr) [200-2750]';
                        'Nb(L) [2750-6000]';
                        'Ni(cr) [200-631]';
                        'Ni(cr) [631-1728]';
                        'Ni(L) [1728-6000]';
                        'P(cr) [195.4-317.3]';
                        'P(L) [317.3-6000]';
                        'Pb(cr) [200-600.65]';
                        'Pb(L) [600.65-3600]';
                        'Rb(cr) [100-312.47]';
                        'Rb(L) [312.47-2100]';
                        'S(a) [200-368.3]';
                        'S(b) [368.3-388.36]';
                        'S(L) [388.36-6000]';
                        'Sc(a) [100-1609]';
                        'Sc(b) [1609-1814]';
                        'Sc(L) [1814-6000]';
                        'Si(cr) [200-1690]';
                        'Si(L) [1690-6000]';
                        'Sn(cr) [200-505.118]';
                        'Sn(L) [505.118-4700]';
                        'Sr(a) [100-820]';
                        'Sr(b) [820-1041]';
                        'Sr(L) [1041-6000]';
                        'Ta(cr) [200-3258]';
                        'Ta(L) [3258-6000]';
                        'Th(a) [200-1650]';
                        'Th(b) [1650-2023]';
                        'Th(L) [2023-6000]';
                        'Ti(a) [200-1156]';
                        'Ti(b) [1156-1944]';
                        'Ti(L) [1944-6000]';
                        'U(a) [200-942]';
                        'U(b) [942-1049]';
                        'U(c) [1049-1408]';
                        'U(L) [1408-4000]';
                        'V(cr) [200-2190]';
                        'V(L) [2190-6000]';
                        'W(cr) [200-3680]';
                        'W(L) [3680-6000]';
                        'Zn(cr) [200-692.73]';
                        'Zn(L) [692.73-6000]';
                        'Zr(a) [200-1135]';
                        'Zr(b) [1135-2125]';
                        'Zr(L) [2125-6000]'};
        end

        function [cp0, hf, h0, ef, s0, g0] = getSpeciesProperties(obj, DB, species, T, MassOrMolar, echo)
            % Calculates the thermodynamic properties of any species included in the NASA database
            %
            % Args:
            %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %     species (char): Chemical species
            %     T (float): Temperature [K]
            %     MassOrMolar (char): Label indicating mass [kg] or molar [mol] units
            %     echo (float): 0 or 1 indicating species not found
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
        
            if nargin < 5, echo = 0; end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Sample application
            %
            % >> [formula, mm, cP0, hf0, h0, ef0, s0, g0] = getSpeciesProperties(obj, DB, 'CO', 1000, 'molar')
            % -------------------------------
            % Possible phases of this species
            % - CO
            % -------------------------------
            %
            % formula =
            %     'C   1.00O   1.00    0.00    0.00    0.00'
            % mm  = 28.0101
            % cP0 = 33.1788
            % hf0 = -1.1054e+05
            % h0  = -8.8848e+04
            % ef0 = -1.1177e+05
            % s0  = 234.5409
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % Change lowercase 'l' to uppercase 'L'
            species = replace(species, 'Al', 'AL');
            species = replace(species, 'Cl', 'CL');
            species = replace(species, 'Tl', 'TL');
            species = replace(species, 'Fl', 'FL');
        
            % Store species name with parenthesis (i.e., the name appearing in NASA's
            % document tables)
            Species_with_parenthesis = species;
        
            % Substitute opening and closing parenthesis and other reserved characters
            % by 'b' in order to format the species name as in the thermo.inp
            % electronic database
        
            FLAG_MILLENIUM = false;
        
            if contains(species, '_M')
                species = strrep(species, '_M', '');
                FLAG_MILLENIUM = true;
            end
        
            name = species;
        
            if name(end) == '+'
                name = [name(1:end - 1) 'plus'];
            elseif name(end) == '-'
                name = [name(1:end - 1) 'minus'];
            end
        
            ind = regexp(name, '[()]');
            name(ind) = 'b';
            ind = regexp(name, '\W');
            name(ind) = '_';
        
            if regexp(name(1), '[0-9]')
                name = ['num_' name];
            end
        
            if FLAG_MILLENIUM
                name = strcat(name, '_M');
            end
        
            species = name;
        
            % If the given species does not exist in strDB, abort the program
            if ~isfield(DB, species)
        
                if echo == 1
                    disp('-------------------------------')
                    disp(['Species ''', species, ''' does not exist as a field in strDB structure'])
                    disp('Program aborted!')
                    disp('-------------------------------')
                end
        
                formula = [];
                W = [];
                cp0 = [];
                hf = [];
                h0 = [];
                ef = [];
                s0 = [];
                g0 = [];
                return
            end
        
            % Detect the position of the phase specifier
            n_open_parenthesis = detect_location_of_phase_specifier(Species_with_parenthesis);
        
            % If it does exist, look for other possible states of aggregation of the
            % same species in DB
            if echo == 1
                % Look for other possible states of aggregation of the same species in
                % strDB
                names = fieldnames(DB);
                any_other_phase = strfind(names, species(1:n_open_parenthesis - 1));
                any_other_phase_index = find(~cellfun(@isempty, any_other_phase));
        
                if ~isempty(any_other_phase_index)
                    disp('-------------------------------')
                    disp('Possible phases of this species')
        
                    for i = 1:length(any_other_phase_index)
                        name_other_phase_with_parenthesis = name_with_parenthesis(names{any_other_phase_index(i)});
                        n_open_parenthesis_other_phase = detect_location_of_phase_specifier(name_other_phase_with_parenthesis);
        
                        if strcmp(name_other_phase_with_parenthesis(1:n_open_parenthesis_other_phase - 1), species(1:n_open_parenthesis - 1))
                            name_other_phase = name_other_phase_with_parenthesis;
                            name_other_phase(name_other_phase(:) == '(') = 'b';
                            name_other_phase(name_other_phase(:) == ')') = 'b';
        
                            if ~iscell(DB.(name_other_phase).tRange)
                                disp(['> ', name_other_phase_with_parenthesis, ' [', num2str(DB.(name_other_phase).tRange(1)), ' K]'])
                            else
                                disp(['> ', name_other_phase_with_parenthesis, ' [', num2str(DB.(name_other_phase).tRange{1}(1)), ' - ', num2str(DB.(name_other_phase).tRange{DB.(name_other_phase).Tintervals}(2)), ' K]'])
                            end
        
                        end
        
                    end
        
                    disp('-------------------------------')
                end
        
            end
        
            % If it does exist, read the corresponding field and store it in the
            % following variables
        
            Tintervals = DB.(species).Tintervals;
            formula = DB.(species).formula;
            phase = DB.(species).phase;
            W = DB.(species).W;
            hf = DB.(species).hf;
            Trange = DB.(species).Trange;
            Texponents = DB.(species).Texponents;
        
            % Set elements and elementsTemperatureReference lists
            [elements, ~] = set_elements(); % Sets elements list
            elementMatrix = set_element_matrix(formula, elements); % Sets elementMatrix
            elementsTemperatureReference = obj.setElementsTemperatureReference(); % Sets elementsTemperatureReference list
        
            % In order to compute the internal energy of formation from the enthalpy of
            % formation of a given species, we must determine the change in moles of
            % gases during the formation reaction of a mole of that species starting
            % from the elements in their reference state. The only elements that are
            % stable as diatomic gases are elements 1 (H), 8 (N), 9 (O), 10 (F), and 18
            % (Cl). The remaining elements that are stable as (monoatomic) gases are
            % the noble gases He (3), Ne (11), Ar (19), Kr (37), Xe (55), and Rn (87),
            % which do not form any compound.
            Delta_n = compute_change_moles_gas_reaction(elementMatrix, phase);
        
            R0 = 8.3144598; % [J/(K-mol)]. Universal gas constant
            % Check if there is at least one temperature interval and, in that case,
            % check that the specified temperature is within limits. If it is not, then
            % abort, otherwise keep on running
            if Tintervals > 0
                a = DB.(species).a;
                b = DB.(species).b;
        
                Tref = obj.temperatureReference;
        
                if (T < Trange{1}(1)) || (T > Trange{Tintervals}(2)) && (echo == 1)
                    disp(['T - out of range [', num2str(Trange{1}(1)), ' - ', num2str(Trange{Tintervals}(2)), ' K] for ', name_with_parenthesis(species)])
                    cp0 = [];
                    h0 = [];
                    ef = hf - Delta_n * R0 * Tref;
                    s0 = [];
                    g0 = [];
                    return
                end
        
                % Get temperature interval
                Tinterval = obj.get_interval(species, T, DB);
                % Compute the thermochemical data at the specified temperature using
                % the polynomial coefficients in the selected temperature interval. All
                % magnitudes are computed in a per mole basis
                cp0 = R0 * sum(a{Tinterval} .* T.^Texponents{Tinterval});
                h0 = R0 * T * (sum(a{Tinterval} .* T.^Texponents{Tinterval} .* [-1 log(T) 1 1/2 1/3 1/4 1/5 0]) + b{Tinterval}(1) / T);
                ef = hf - Delta_n * R0 * Tref;
                s0 = R0 * (sum(a{Tinterval} .* T.^Texponents{Tinterval} .* [-1/2 -1 log(T) 1 1/2 1/3 1/4 0]) + b{Tinterval}(2));
        
                % Compute the standar gibbs free energy of formation at the specified
                % temperature. This enforces us to consider explicitely the formation
                % reaction from the elements in their reference states at room
                % temperature, unless the species is precisely an element in its
                % reference state, in which case the standard gibbs free energy of
                % formation is identically zero.
        
                % disp(['Species = ',Species])
                % name_with_parenthesis(Species)
        
                [iRE, REname] = isRefElm(elementsTemperatureReference, species(1:n_open_parenthesis - 1), T);
        
                if (~iRE)
        
                    if echo == 1
                        disp([species, ' is not Ref-Elm.'])
                    end
        
                    g0 = h0 - T .* s0;
                else
        
                    if echo == 1
                        disp([REname, ' is Ref-Elm.'])
                    end
        
                    g0 = 0;
                end
        
                if strcmpi(MassOrMolar, 'mass')
                    cp0 = molar2mass(cp0, W);
                    hf = molar2mass(hf, W);
                    ef = molar2mass(ef, W);
                    h0 = molar2mass(h0, W);
                    s0 = molar2mass(s0, W);
        
                    if phase == 0
                        g0 = molar2mass(g0, W);
                    else
                        g0 = [];
                    end
        
                end
        
            % If the species is only a reactant determine it's reference temperature
            % Tref. For noncryogenic reactants, assigned enthalpies are given at 298.15
            % K. For cryogenic liquids, assigned enthalpies are given at their boiling
            % points instead of 298.15 K
            else
        
                if T ~= Trange(1)
                    disp(['T - out of range for ', name_with_parenthesis(species), ' [', num2str(Trange(1)), ' K]'])
                    cp0 = [];
                    h0 = [];
                    ef = hf - Delta_n * R0 * Trange(1);
                    s0 = [];
                    g0 = [];
        
                    if strcmpi(MassOrMolar, 'mass')
                        hf = molar2mass(hf, W);
                        ef = molar2mass(ef, W);
                    end
        
                    return
                end
        
                Tref = Trange(1);
        
                cp0 = 0;
                h0 = hf;
                ef = hf - Delta_n * R0 * Tref;
                s0 = 0;
                g0 = hf;
        
                if strcmpi(MassOrMolar, 'mass')
                    hf = molar2mass(hf, W);
                    ef = molar2mass(ef, W);
                end
        
            end

            % SUB-PASS FUNCTIONS
            function value = molar2mass(value, W)
                % Change molar units [mol] to mass units [kg]
                value = value / W * 1e3;
            end

            function Delta_n = compute_change_moles_gas_reaction(elementMatrix, phase)
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
                %     elementMatrix (float): Element matrix of the species
                %     phase (float): 0 or 1 indicating gas or condensed species
                %
                % Returns:
                %     Delta_n (float): Change in moles of gases during the formation reaction of a mole of that species starting from the elements in their reference state
                %
                % Example:
                %     Delta_n = compute_change_moles_gas_reaction(element_matrix, phase)
            
                Delta_n_per_mole = sum(elementMatrix(1,:) == [1, 8, 9, 10, 18]') / 2 ... 
                                 + sum(elementMatrix(1,:) == [3, 11, 19, 37, 55, 87]');
                Delta_n = 1 - phase - dot(Delta_n_per_mole, elementMatrix(2,:));
            end

            function [FLAG_RE, REname] = isRefElm(reference_elements, species, T)
                % Check if the given species is a reference element
                %
                % Args:
                %     reference_elements (cell): List of reference elements with temperature intervals [K]
                %     species (char): Chemical species
                %     T (float): Temperature
                %
                % Returns:
                %     name (char): Full name of the given species
                %
                % Example:
                %     [FLAG_RE, REname] = isRefElm(reference_elements, 'O', 1000)
            
                % Change lowercase 'l' to uppercase 'L' for Al, Cl, Tl, and Fl
                species(strfind(species, 'Al') + 1) = 'L';
                species(strfind(species, 'Cl') + 1) = 'L';
                species(strfind(species, 'Tl') + 1) = 'L';
                species(strfind(species, 'Fl') + 1) = 'L';
            
                FLAG_RE = false;
                REname = [];
            
                % Look for entries in the Reference_form_of_elements_with_T_intervals list
                % that partially match with the desired species and then check each one
                % sucessivelly
            
                j = find(contains(reference_elements, species));
            
                for i = 1:length(j)
                    % disp(num2str(i))
                    TentativeRefElm = reference_elements{j(i)};
                    % Detect temperature interval
                    n1 = strfind(TentativeRefElm, '[');
                    n2 = strfind(TentativeRefElm, '-');
                    n3 = strfind(TentativeRefElm, ']');
                    T1 = sscanf(TentativeRefElm(n1 + 1:n2 - 1), '%f');
                    T2 = sscanf(TentativeRefElm(n2 + 1:n3 - 1), '%f');
            
                    if (T < T1) && (T > T2)
                        continue
                    end

                    % Detect location of open parenthesis
                    n_open_parenthesis = strfind(TentativeRefElm(1:n1 - 2), '(');
                    % Detect location of '2'
                    n_two = strfind(TentativeRefElm(1:n1 - 2), '2');
                    % If thera are no '2's or parenthesis, the Species is essentially a noble gas
                    if isempty(n_open_parenthesis) && isempty(n_two)
                        % 1
                        if strcmp(TentativeRefElm(1:n1 - 2), species)
                            FLAG_RE = true;
                            REname = TentativeRefElm(1:n1 - 2);
                        end
        
                    end
        
                    % If there are '2's the species may be in the reference state or
                    % not (e.g. O2 is, but O is not)
                    if ~isempty(n_two)
                        % 2
                        if strcmp(TentativeRefElm(1:n_two), species)
                            FLAG_RE = true;
                            REname = TentativeRefElm(1:n_two);
                        end
        
                        if strcmp(TentativeRefElm(1:n_two - 1), species)
                            REname = TentativeRefElm(1:n_two);
                        end
        
                    end
        
                    % If there are opening parenthesis, the species is in condensed phase
                    if ~isempty(n_open_parenthesis)
                        % 3
                        if strcmp(TentativeRefElm(1:n_open_parenthesis - 1), species)
                            FLAG_RE = true;
                            REname = TentativeRefElm(1:n1 - 2);
                        end
        
                    end
            
                end
            
            end

        end

        function Tinterval = get_interval(species, T, DB)
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

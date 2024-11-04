classdef (Abstract) Database < handle
    % The :mat:class:`Database` abstract class contains the common methods between database objects.
    % This class is used as a base class for the :mat:class:`NasaDatabase` and :mat:class:`BurcatDatabase` classes.
    %
    % See also: :mat:func:`NasaDatabase`, :mat:func:`BurcatDatabase`

    properties
        name                   % Database name
        species                % Struct with Species objects
        filename               % Database filename
        interpolationMethod    % Interpolation method for the griddedInterpolant objects
        extrapolationMethod    % Extrapolation method for the griddedInterpolant objects
        pointsTemperature      % Number of points to use in the interpolation of the thermodynamic data
        temperatureReference   % Default temperature of reference
        thermoFile             % Name with path for the thermodynamic database file (raw data)
        units                  % Units thermodynamic data [molar or mass]
    end

    properties (Hidden)
        id                     % Database ID
        FLAG_BENCHMARK = false % Flag to reload the database despite being in cache for benchamarking purposes
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
            addParameter(ip, 'interpolationMethod', defaultInterpolationMethod, @ischar);
            addParameter(ip, 'extrapolationMethod', defaultExtrapolationMethod, @ischar);
            addParameter(ip, 'units', defaultUnits, @ischar);
            addParameter(ip, 'pointsTemperature', defaultPointsTemperature, @isnumeric);
            addParameter(ip, 'temperatureReference', defaultTemperatureReference, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
            addParameter(ip, 'thermoFile', defaultThermoFile);
            addParameter(ip, 'FLAG_BENCHMARK', obj.FLAG_BENCHMARK);
            parse(ip, varargin{:});
            
            % Set properties
            obj.name = ip.Results.name;
            obj.species = ip.Results.species;
            obj.filename = ip.Results.filename;
            obj.interpolationMethod = ip.Results.interpolationMethod;
            obj.extrapolationMethod = ip.Results.extrapolationMethod;
            obj.units = ip.Results.units;
            obj.pointsTemperature = ip.Results.pointsTemperature;
            obj.temperatureReference = ip.Results.temperatureReference;
            obj.thermoFile = ip.Results.thermoFile;

            % Set ID
            setID(obj);

            % Check if database is in cached and the id matches
            if ip.Results.FLAG_BENCHMARK
                obj = obj.load();
                return
            end

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
            % Load database from file
            %
            % Args:
            %     obj (Database): Database object
            %
            % Optional Args:
            %     filename (char): Filename of the database
            %
            % Returns:
            %     obj (Database): Database object

            if nargin > 1
                obj.filename = varargin{1};
            end

            % Load database
            if exist(obj.filename, 'file')
                fprintf('%s database with thermo loaded from the main path ... ', obj.name);
                load(obj.filename, 'DB');
                obj = DB;
            else
                generateDatabase(obj);
            end

            % Status
            fprintf('OK!\n');
        end

        function save(obj)
            % Save database to a *.mat file
            % 
            % Args:
            %     obj (Database): Database object

            DB = obj;
            uisave({'DB'}, 'DB.mat');
        end

        function value = getProperty(obj, listSpecies, property)
            % Gets the vector of the defined property for the given
            % set of species
            %
            % Args:
            %     obj (Database): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
            %     listSpecies (cell): List of species
            %     property (char): Property name
            %     
            %
            % Returns:
            %     value (float): Property vector
            %
            % Example:
            %     value = getProperty(obj, {'H2O', 'CO2'}, 'hf')
        
            for i = length(listSpecies):-1:1
                value(1, i) = obj.species.(listSpecies{i}).(property);
            end
            
        end

    end

    methods (Access = protected)

        function obj = setID(obj)
            % Concatenate input arguments to create a unique identifier string
            value =  [obj.name, num2str(obj.species), obj.filename, ...
                      obj.interpolationMethod, obj.extrapolationMethod, obj.units, ...
                      num2str(obj.pointsTemperature), num2str(obj.temperatureReference),...
                      obj.thermoFile];

            obj.id = combustiontoolbox.utils.generateID(value);
        end

    end

    methods (Abstract)
        generateDatabase(obj)
        getSpeciesThermo(obj)
    end

    methods (Access = public)
        
        function addSpecies(obj, speciesName, DB_master)
            % Add species to the database

            speciesName = obj.fullname2name(speciesName);
    
            % Initialization
            species = DB_master.(speciesName);
    
            % Get data
            Tintervals = species.Tintervals;
            Trange = species.Trange;
    
            if Tintervals == 0
                % Handle species with no temperature intervals
                species = obj.computeConstantTemperatureSpecies(species, speciesName, Trange, DB_master);
            else
                % Handle species with temperature intervals
                species = obj.computeVariableTemperatureSpecies(species, speciesName, Trange, Tintervals, DB_master);
            end
    
            % Store the species data in the species struct property
            obj.species.(speciesName) = species;
        end
        
    end

    methods (Access = private)

        function species = computeConstantTemperatureSpecies(obj, species, speciesName, Trange, DB_master)
            Tref = Trange(1);
        
            % Get thermodynamic data at reference temperature
            [Cp0, Hf0, H0, Ef0, S0, DfG0] = obj.getSpeciesThermo(DB_master, speciesName, Tref, obj.units);
        
            % Store the thermodynamic data
            species.hf = Hf0;
            species.ef = Ef0;
            species.Tref = Tref;
            species.T = Tref;
        
            % Generate interpolation curves (constant values)
            species.cpcurve = griddedInterpolant([Tref, Tref + 1], [Cp0, Cp0], 'linear', 'linear');
            species.h0curve = griddedInterpolant([Tref, Tref + 1], [H0, H0], 'linear', 'linear');
            species.s0curve = griddedInterpolant([Tref, Tref + 1], [S0, S0], 'linear', 'linear');
            species.g0curve = griddedInterpolant([Tref, Tref + 1], [DfG0, DfG0], 'linear', 'linear');
        end
        
        function species = computeVariableTemperatureSpecies(obj, species, speciesName, Trange, Tintervals, DB_master)
            Tmin = Trange{1}(1);
            Tmax = Trange{Tintervals}(2);
            T_vector = linspace(Tmin, Tmax, obj.pointsTemperature);
            
            % Store the thermodynamic data
            species.T = T_vector;
            [~, Hf0, ~, Ef0, ~, ~] = obj.getSpeciesThermo(DB_master, speciesName, species.Tref, obj.units);
            species.hf = Hf0;
            species.ef = Ef0;
            
            % Get thermodynamic data over the temperature range
            [Cp0_vector, ~, H0_vector, ~, S0_vector, ~] = obj.getSpeciesThermo(DB_master, speciesName, T_vector, obj.units);
            DfG0_vector = H0_vector - T_vector .* S0_vector;
        
            % Generate interpolation curves
            species.cpcurve = griddedInterpolant(T_vector, Cp0_vector, obj.interpolationMethod, obj.extrapolationMethod);
            species.h0curve = griddedInterpolant(T_vector, H0_vector, obj.interpolationMethod, obj.extrapolationMethod);
            species.s0curve = griddedInterpolant(T_vector, S0_vector, obj.interpolationMethod, obj.extrapolationMethod);
            species.g0curve = griddedInterpolant(T_vector, DfG0_vector, obj.interpolationMethod, obj.extrapolationMethod);
        
            % Store additional species data
            species.Tintervals = species.Tintervals;
            species.Trange = species.Trange;
            species.Texponents = species.Texponents;
            species.a = species.a;
            species.b = species.b;
        end

    end

end

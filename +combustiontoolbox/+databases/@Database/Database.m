classdef (Abstract) Database < handle
    % The :mat:class:`Database` abstract class contains the common methods between database objects.
    % This class is used as a base class for the :mat:class:`NasaDatabase` and :mat:class:`BurcatDatabase` classes.
    %
    % See also: :mat:func:`NasaDatabase`, :mat:func:`BurcatDatabase`

    properties
        name
        species
        filename
        interpolationMethod
        extrapolationMethod 
        units
        pointsTemperature
        temperatureReference
        thermoFile
    end

    properties (Hidden)
        id
        FLAG_BENCHMARK = false
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

            % Generate ID
            obj = obj.generate_id();

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

        function obj = generate_id(obj)
            % Concatenate input arguments to create a unique identifier string
            value =  [obj.name, num2str(obj.species), obj.filename, ...
                      obj.interpolationMethod, obj.extrapolationMethod, obj.units, ...
                      num2str(obj.pointsTemperature), num2str(obj.temperatureReference),...
                      obj.thermoFile];

            obj.id = combustiontoolbox.utils.generate_id(value);
        end

    end
    
end

classdef (Abstract) Database < handle

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
                DB = obj.generateDatabase(DB_master);
                obj.species = DB;
            end

            % Status
            fprintf('OK!\n');
        end

        function save(obj)
            % Save database
            DB = obj;
            DB_master = obj.species;
            uisave({'DB'}, 'DB.mat');
            uisave({'DB_master'}, 'DB_master.mat');
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
    
end

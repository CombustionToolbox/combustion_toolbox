classdef BenchmarkMetadata < handle
    % The :mat:func:`BenchmarkMetadata` class is used to store and manage metadata
    % information for benchmark tests in the Combustion Toolbox.
    %
    % The :mat:func:`BenchmarkMetadata` object can be initialized as follows: ::
    %
    %       metadata = BenchmarkMetadata(solver, mixtureArray, 'filename', 'run_validation_TP_CEA_1', 'testname', 'C6H6_air1_TP')
    %
    % Here the parameters specify various metadata attributes such as the module name,
    % problem type, filename, test name, number of cases, number of species, and tolerance.
    %
    % See also: :mat:func:`Benchmark`

    properties
        module       % Name of the module being benchmarked
        solver       % Solver used in the benchmark
        problemType  % Type of problem (e.g., 'TP', 'TV', 'HP', etc.)
        filename     % Name of the file containing benchmark data
        numCases     % Number of cases in the benchmark
        numSpecies   % Number of species involved
        tolerance    % Tolerance level for the benchmark
        avgTime      % Average execution time for the benchmark
    end

    properties (Access = private)
        mixtureArray % Array of mixture objects used in the benchmark
    end

    methods
        function obj = BenchmarkMetadata(solver, mixtureArray, filename, varargin)
            % Constructor

            % Parse input arguments
            p = inputParser;
            addRequired(p, 'solver', @(x) isobject(x));
            addRequired(p, 'mixtureArray', @(x) isobject(x));
            addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
            parse(p, solver, mixtureArray, filename, varargin{:});

            % Set properties
            obj.solver = p.Results.solver;
            obj.mixtureArray = p.Results.mixtureArray;
            obj.filename = p.Results.filename;
            
            % Get metadata from solver and mixtureArray
            obj = obj.getMetadata();
        end

        function obj = setTime(obj, time)
            % Set average execution time
            obj.avgTime = time;
        end

        function print(obj)
            % Display metadata information
            fprintf('Benchmark Metadata:\n');
            fprintf('  Module:            %s\n', obj.module);
            fprintf('  Problem Type:      %s\n', obj.problemType);
            fprintf('  Filename:          %s\n', obj.filename);
            fprintf('  Test Name:         %s\n', obj.testname);
            fprintf('  Number of Cases:   %d\n', obj.numCases);
            fprintf('  Number of Species: %d\n', obj.numSpecies);
            fprintf('  Tolerance:         %s\n', num2str(obj.tolerance));
            fprintf('  Average Time:      %s seconds\n', num2str(obj.time));
        end

        function T = asTable(obj)
            % Export metadata as a single-row MATLAB table

            % Get clean class name (no package)
            solverClass = string(metaclass(obj.solver).Name);
            if contains(solverClass, '.')
                parts = strsplit(solverClass, '.');
                solverClass = parts{end};
            end

            tolStr = arrayfun(@(x) sprintf('%.2e', x), obj.tolerance, 'UniformOutput', false);

            T = table( ...
                categorical(string(obj.module)), ...
                categorical(string(solverClass)), ...
                categorical(string(obj.problemType)), ...
                categorical(string(obj.filename)), ...
                uint32(obj.numCases), ...
                uint32(obj.numSpecies), ...
                categorical(tolStr), ...
                obj.avgTime, ...
                'VariableNames', {'Module','Solver','Problem','Filename','Cases','Species','Tolerance','AvgTime'} ...
            );
        end

    end

    methods (Access = private)
        function obj = getMetadata(obj)
            % Get metadata

            % Get module
            switch class(obj.solver)
                case {'combustiontoolbox.equilibrium.EquilibriumSolver'}
                    obj.module = 'CT-EQUIL';
                case {'combustiontoolbox.shockdetonation.ShockSolver', 'combustiontoolbox.shockdetonation.DetonationSolver', 'combustiontoolbox.shockdetonation.JumpConditionsSolver'}
                    obj.module = 'CT-SD';
                case {'combustiontoolbox.rocket.RocketSolver'}
                    obj.module = 'CT-ROCKET';
                case {'combustiontoolbox.turbulence.HelmholtzSolver'}
                    obj.module = 'CT-TURBULENCE';
                otherwise
                    obj.module = 'Unknown';
            end

            % Get problem type
            obj.problemType = obj.solver.problemType;

            % Get number of cases
            obj.numCases = length(obj.mixtureArray);

            % Get number of species
            obj.numSpecies = length(obj.mixtureArray(1).chemicalSystem.listSpecies);

            % Get tolerance
            if isa(obj.solver, 'combustiontoolbox.equilibrium.EquilibriumSolver') && any(strcmp(obj.problemType, {'TP','TV'}))
                obj.tolerance = obj.solver.tolGibbs;
            elseif isprop(obj.solver, 'tol0')
                obj.tolerance = obj.solver.tol0;
            else
                obj.tolerance = NaN;
            end

        end

    end

end
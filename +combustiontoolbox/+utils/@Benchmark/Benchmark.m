classdef Benchmark < handle
    % The :mat:func:`Benchmark` class is used to perform a set of benchmark tests using the 
    % Combustion Toolbox code. It measures and reports the average computational time of 
    % selected functions over multiple iterations.
    % 
    % The :mat:func:`Benchmark` object can be initialized as follows: ::
    %
    %       bench = Benchmark(tests, numIterations)
    %
    % Here ``tests`` is a cell array of function handles representing the benchmarked functions,
    % and ``numIterations`` specifies the number of times each function will be executed to 
    % compute the average execution time.
    %
    % Example usage:
    %     % Define benchmark tests as function handles
    %     tests = {@run_validation_TP_CEA_6, @another_test_function};
    %
    %     % Create the Benchmark object with tests and number of iterations
    %     bench = Benchmark(tests, 20);
    %
    %     % Run the tests
    %     bench = bench.run();
    %
    %     % Generate and display the report
    %     bench.report();
    %
    % See also: :mat:func:`combustiontoolbox.utils.timeFunction`, :mat:func:`combustiontoolbox.utils.extensions.cpuinfo`


    properties
        tests         % Cell array of function handles
        numIterations % Number of iterations for each test
        results       % Struct to store results
        system        % Struct storing CPU and system information
    end

    properties (Dependent)
        numTest       % Number of tests 
    end
    
    methods

        function obj = Benchmark(varargin)
            % Constructor
            defaultTests = {@run_validation_TP_CEA_6, @run_validation_DET_CEA_1, @run_validation_SHOCK_POLAR_SDToolbox_1, @run_validation_ROCKET_CEA_23};
            defaultNumIterations = 10;

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'tests', defaultTests, @(x) iscell(x));
            addParameter(p, 'numIterations', defaultNumIterations, @(x) isnumeric(x) && x > 0);
            parse(p, varargin{:});

            % Set properties
            obj.tests = p.Results.tests;
            obj.numIterations = p.Results.numIterations;
            obj.results = struct();
            obj.system = combustiontoolbox.utils.extensions.cpuinfo();
        end

        function value = get.numTest(obj)
            % Get number of tests
            value = length(obj.tests);
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the EquilibriumSolver object
            %
            % Args:
            %     obj (Benchmark): Benchmark object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (Benchmark): Benchmark object with updated properties
            %
            % Examples:
            %     * set(Benchmark(), 'tests', {@run_validation_TP_CEA_5, @run_validation_TP_CEA_6});
            %     * set(Benchmark(), 'numIterations', 20);
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end
        
        function obj = run(obj)
            % The :mat:func:`run` method executes the benchmark tests and records the average 
            % execution time for each test function.
            %
            % Args:
            %     obj (Benchmark): Benchmark object containing the tests to be executed.
            %
            % Returns:
            %     obj (Benchmark): Benchmark object updated with the execution times of each test.
            %
            % Examples:
            %     * bench = bench.run();

            for i = 1:obj.numTest
                testFunc = obj.tests{i};
                testName = func2str(testFunc);
                
                fprintf('Running test: %s...\n', testName);
                
                avgTime = combustiontoolbox.utils.timeFunction(testFunc, obj.numIterations);
                obj.results.(testName) = avgTime;
                
                fprintf('%-30s | Average Time = %.6f seconds\n', testName, avgTime);
            end

        end
        
        function report(obj)
            % The :mat:func:`report` method displays a formatted benchmark report, showing the 
            % average execution times for each test.
            %
            % Args:
            %     obj (Benchmark): Benchmark object containing recorded execution times.
            %
            % Examples:
            %     * bench.report();
            %
            % See also: :mat:func:`Benchmark`, :mat:func:`runTests`

            if isempty(fieldnames(obj.results))
                warning('No benchmark results available. Run the benchmarks first.');
                return;
            end
            
            % Display CPU/system summary
            fprintf('************************************************\n');
            fprintf('------------------------------------------------\n');
            fprintf('\nSystem Information:\n');
            fprintf('------------------------------------------------\n');
            fprintf('CPU Name:              %s\n', obj.system.CPUName);
            fprintf('Clock Speed:           %s\n', obj.system.Clock);
            fprintf('CPU Cache (L2):        %d bytes\n', obj.system.Cache);
            fprintf('CPU Cores:             %d\n', obj.system.TotalCores);
            fprintf('Operating System:      %s %s\n', obj.system.OSType, obj.system.OSVersion);
            fprintf('Hostname:              %s\n', obj.system.Hostname);
            fprintf('------------------------------------------------\n');

            % Display benchmark results
            fprintf('\nBenchmark Report:\n');
            fprintf('------------------------------------------------\n');
            testNames = fieldnames(obj.results);

            for i = 1:obj.numTest
                fprintf('%-30s | %.6f s\n', testNames{i}, obj.results.(testNames{i}));
            end

            fprintf('------------------------------------------------\n');
            fprintf('************************************************\n');
        end
    end
end

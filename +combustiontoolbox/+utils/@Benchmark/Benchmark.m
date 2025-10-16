classdef Benchmark < handle
    % The :mat:func:`Benchmark` class is used to perform a set of benchmark
    % tests using the Combustion Toolbox code. It measures and reports the 
    % average computational time of selected functions over multiple 
    % iterations.
    % 
    % The :mat:func:`Benchmark` object can be initialized as follows: ::
    %
    %       bench = Benchmark(tests, numIterations)
    %
    % Here `tests` is a cell array of function handles representing the 
    % benchmarked functions, and `numIterations` specifies the number of 
    % times each function will be executed to compute the average execution
    % time.
    %
    % Example usage:
    %     % Define benchmark tests as function handles
    %     tests = {@run_validation_TP_CEA_6, @another_test_function};
    %
    %     % Create the Benchmark object with tests and number of iterations
    %     bench = Benchmark(tests, 'numIterations', 20);
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
        metadata      % BenchmarkMetadata array storing metadata for each test
        system        % Struct storing CPU and system information
    end

    properties (Dependent)
        numTest       % Number of tests 
    end
    
    methods

        function obj = Benchmark(varargin)
            % Constructor
            defaultTests = obj.getDefaultTests();
            defaultNumIterations = 10;

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'tests', defaultTests, @(x) iscell(x));
            addParameter(p, 'numIterations', defaultNumIterations, @(x) isnumeric(x) && x > 0);
            parse(p, varargin{:});

            % Set properties
            obj.tests = p.Results.tests;
            obj.numIterations = p.Results.numIterations;
            obj.metadata = combustiontoolbox.utils.BenchmarkMetadata.empty();
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

            % Import packages
            import combustiontoolbox.utils.timeFunction

            for i = 1:obj.numTest
                testFunc = obj.tests{i};
                testName = func2str(testFunc);
                
                fprintf('Running test: %s...\n', testName);

                [avgTime, ~, metadata] = timeFunction(testFunc, obj.numIterations, 'FLAG_BENCHMARK', true); %#ok<PROP>
                
                obj.metadata(end + 1) = metadata; %#ok<PROP>
                obj.metadata(end) = obj.metadata(end).setTime(avgTime);
                
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
            
            % Import packages
            import combustiontoolbox.utils.prettyHeader
            
            if isempty(fieldnames(obj.metadata))
                warning('No benchmark results available. Run the benchmarks first.');
                return;
            end

            % Get Combustion Toolbox version
            vCT = erase(combustiontoolbox.common.Constants.release, 'v');

            % Get MATLAB version
            vMATLAB = matlabRelease().Release;

            % Handle reformat to categorical
            toCat = @(x) categorical(string(x));

            % Get CPU/system summary
            table1 = table( ...
                toCat(vCT), toCat(vMATLAB), toCat(obj.system.OSVersion), toCat(obj.system.CPUName),...
                toCat(obj.system.Clock), toCat(obj.system.Cache), toCat(obj.system.TotalCores), ...
                'VariableNames', {'CT version', 'MATLAB version', 'OS version', 'CPU name', 'Clock speed', 'Cache L2', 'Cores'} ...
            );
        
            % Collect tables for each BenchmarkMetadata
            tables = arrayfun(@(m) m.asTable(), obj.metadata, 'UniformOutput', false);

            % Get extended benchmark table
            table2 = vertcat(tables{:});

            % Get summary benchmark table
            table3 = groupsummary(table2, {'Module', 'Problem'}, {'sum', 'mean'}, ...
                                  {'Cases', 'Species', 'AvgTime'});
        
            % Print
            fprintf(sprintf('\n%s\n\n', prettyHeader('BENCHMARK REPORT', 'unicode')));

            fprintf(sprintf('\n%s\n\n', prettyHeader('System Information', 'unicode')));
            disp(table1);

            fprintf(sprintf('%s\n\n', prettyHeader('Benchmarking', 'unicode')));
            disp(table2);
        
            fprintf(sprintf('%s\n\n', prettyHeader('Summary Benchmarking', 'unicode')));
            disp(table3);
        end

    end

    methods (Static)
        
        function defaultTests = getDefaultTests()
            % Get default tests
            defaultTests = {...
                @run_validation_TP_TEA_1,...
                @run_validation_TP_TEA_2,...
                @run_validation_TP_TEA_3,...
                @run_validation_TP_TEA_4,...
                @run_validation_TP_CEA_1,...
                @run_validation_TP_CEA_2,...
                @run_validation_TP_CEA_3,...
                @run_validation_TP_CEA_4,...
                @run_validation_TP_CEA_5,...
                @run_validation_TP_CEA_6,...
                @run_validation_TV_CEA_1,...
                @run_validation_TV_CEA_2,...
                @run_validation_HP_CEA_1,...
                @run_validation_HP_CEA_2,...
                @run_validation_HP_CEA_3,...
                @run_validation_HP_CEA_4,...
                @run_validation_EV_CEA_1,...
                @run_validation_SP_CEA_1,...
                @run_validation_SV_CEA_1,...
                @run_validation_SHOCK_IONIZATION_CEA_1,...
                @run_validation_SHOCK_IONIZATION_CEA_2,...
                @run_validation_SHOCK_R_IONIZATION_CEA_1,...
                @run_validation_SHOCK_R_IONIZATION_CEA_2,...
                @run_validation_SHOCK_POLAR_SDToolbox_1,...
                @run_validation_SHOCK_POLAR_SDToolbox_2,...
                @run_validation_SHOCK_PRANDTL_MEYER_SDToolbox_1, ...
                @run_validation_SHOCK_PRANDTL_MEYER_SDToolbox_2, ...
                @run_validation_DET_CEA_1,...
                @run_validation_DET_CEA_2,...
                @run_validation_DET_CEA_3,...
                @run_validation_DET_CEA_4,...
                @run_validation_DET_OVERDRIVEN_SDToolbox_1,...
                @run_validation_DET_POLAR_SDToolbox_1,...
                @run_validation_ROCKET_CEA_20,...
                @run_validation_ROCKET_CEA_21,...
                @run_validation_ROCKET_CEA_22,...
                @run_validation_ROCKET_CEA_23,...
                @run_validation_ROCKET_CEA_24,...
                @run_validation_ROCKET_CEA_25,...
            };
        end
        
    end

end
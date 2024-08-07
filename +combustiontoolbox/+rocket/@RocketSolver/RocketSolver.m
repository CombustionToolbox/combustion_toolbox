classdef RocketSolver < handle
    % The :mat:func:`RocketSolver` class is used to solve rocket problems. The class provides methods
    % to compute the performance of a rocket using the Infinite Area Chamber (IAC) and Finite Area
    % Chamber (FAC) models.
    %
    % The :mat:func:`RocketSolver` object can be initialized as follows: ::
    %
    %       solver = RocketSolver('problemType', problemType, ...)
    %
    % Here ``problemType`` represents the acronym of the problem to be solved (see below).
    % Additional optional parameters can be provided to customize the solver's behavior.
    %
    % Problem types:
    %     * ``ROCKET_IAC``: Get the performance of a rocket using the Infinite Area Chamber (IAC) model
    %     * ``ROCKET_FAC``: Get the performance of a rocket using the Finite Area Chamber (FAC) model
    %
    % See also: :mat:func:`Mixture`, :mat:func:`EquilibriumSolver`, :mat:func:`solve`, :mat:func:`solveArray`, :mat:func:`report`

    properties
        problemType           % Problem type
        equilibriumSolver     % EquilibriumSolver object
        % * Rocket propellant performance (CT-ROCKET module)
        tol0 = 1e-4;          % Tolerance rocket performance
        itMax = 10;           % Max number of iterations - rocket performance
        % * Flags
        FLAG_SUBSONIC = false % Flag to indicate subsonic Area ratio
        FLAG_RESULTS = true   % Flag to print results
        FLAG_TIME = true      % Flag to print elapsed time
        FLAG_REPORT = false   % Flag to print predefined plots
        % * Miscellaneous
        time                  % Elapsed time [s]
        plotConfig            % PlotConfig object
    end

    methods

        function obj = RocketSolver(varargin)
            % Constructor
            defaultProblemType = 'ROCKET_IAC';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'ROCKET_IAC', 'ROCKET_FAC'})));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addParameter(p, 'FLAG_TCHEM_FROZEN', defaultFLAG_TCHEM_FROZEN, @(x) islogical(x))
            addParameter(p, 'FLAG_FROZEN', defaultFLAG_FROZEN, @(x) islogical(x));
            addParameter(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            addParameter(p, 'tolMoles', defaultEquilibriumSolver.tolMoles, @(x) isnumeric(x) && x > 0);
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;
            
            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
                obj.equilibriumSolver.FLAG_FROZEN = p.Results.FLAG_FROZEN;
                obj.equilibriumSolver.tolMoles = p.Results.tolMoles;
            end

            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
            obj.plotConfig.plotProperties = [obj.plotConfig.plotProperties, {'u', 'I_sp', 'I_vac'}];
            obj.plotConfig.plotPropertiesBasis = [obj.plotConfig.plotPropertiesBasis, {[], [], []}];
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the RocketSolver object
            %
            % Args:
            %     obj (RocketSolver): RocketSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (RocketSolver): RocketSolver object with updated properties
            %
            % Examples:
            %     * set(RocketSolver(), 'tol0', 1e-6);
            %     * set(RocketSolver(), 'problemType', 'ROCKET_FAC');
            
            varargin = [property, value, varargin];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function varargout = solve(obj, mix1, varargin)
            % Solve rocket problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): initial Mixture object
            %
            % Optional Args:
            %     * mix2_inj_guess (Mixture): Initial guess for the injector Mixture object
            %     * mix2_c_guess (Mixture): Initial guess for the chamber Mixture object
            %     * mix3_guess (Mixture): Initial guess for the throat Mixture object
            %     * mix4_guess (Mixture): Initial guess for the exit Mixture object
            %
            % Returns:
            %     varargout (Mixture): Updated Mixture objects depending on the problem type
            %
            % Examples:
            %     * [mix1, mix2_c, mix3] = solve(RocketSolver(), mix1); % Rocket IAC
            %     * [mix1, mix2_c, mix3, mix4] = solve(RocketSolver(), mix1); % Rocket IAC
            %     * [mix1, mix2_c, mix3] = solve(RocketSolver(), mix1, mix2Guess, mix3Guess); % Rocket IAC
            %     * [mix1, mix2_c, mix3, mix4] = solve(RocketSolver(), mix1, mix2Guess, mix3Guess, mix4Guess); % Rocket IAC
            %     * [mix1, mix2_inj, mix2_c, mix3] = solve(RocketSolver(), mix1); % Rocket FAC
            %     * [mix1, mix2_inj, mix2_c, mix3, mix4] = solve(RocketSolver(), mix1); % Rocket FAC
            %     * [mix1, mix2_inj, mix2_c, mix3] = solve(RocketSolver(), mix1, mix2Guess, mix3Guess); % Rocket FAC
            %     * [mix1, mix2_inj, mix2_c, mix3] = solve(RocketSolver(), mix1, mix2Guess, mix3Guess, mix4Guess); % Rocket FAC
            
            switch upper(obj.problemType)
                case 'ROCKET_IAC'
                    % Solve rocket problem with Infinite Area Chamber (IAC) model
                    if nargin > 2
                        [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, varargin{:});
                    else
                        [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1);
                    end

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2_c.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;
                    mix4.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2_c, mix3, mix4);
                    end

                    % Set output
                    varargout = {mix1, mix2_c, mix3, mix4};
                case 'ROCKET_FAC'
                    % Solve rocket problem considering an Finite Area Chamber (FAC) model
                    if nargin > 2
                        [mix1, mix2_inj, mix2_c, mix3, mix4] = rocketFAC(obj, mix1, varargin{:});
                    else
                        [mix1, mix2_inj, mix2_c, mix3, mix4] = rocketFAC(obj, mix1);
                    end
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2_inj.problemType = obj.problemType;
                    mix2_c.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;
                    mix4.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2_inj, mix2_c, mix3, mix4);
                    end

                    % Set output
                    varargout = {mix1, mix2_inj, mix2_c, mix3, mix4};
                otherwise
                    error('Invalid problem type');
            end

        end

        function varargout = solveArray(obj, mix1Array, varargin)
            % Solve a set of rocket problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1Array (Mixture): Array of initial Mixture objects
            %
            % Returns:
            %     varargout (Mixture): Updated arrays of Mixture objects depending on the shock problem type
            %
            % Examples:
            %     * [mix1Array, mix2Array, mix3Array] = solveArray(RocketSolver(), mix1Array); % Rocket IAC
            %     * [mix1Array, mix2Array, mix3Array, mix4Array] = solveArray(RocketSolver(), mix1Array); % Rocket IAC
            %     * [mix1Array, mix2Array, mix3Array, mix4Array] = solveArray(RocketSolver(), mix1Array); % Rocket FAC
            %     * [mix1Array, mix2Array, mix3Array, mix4Array, mix5Array] = solveArray(RocketSolver(), mix1Array); % Rocket FAC

            % Definitions
            n = length(mix1Array);
            problem = obj.problemType;
            
            % Timer
            obj.time = tic;

            % Initialization
            mix2Array = mix1Array;

            % Calculations
            switch upper(problem)
                case {'ROCKET_IAC'}
                    % Initialization
                    mix3Array = mix1Array;
                    mix4Array = mix1Array;
                    
                    % Calculations
                    [mix1Array(n), mix2Array(n), mix3Array(n), mix4Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i), mix3Array(i), mix4Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1), mix3Array(i + 1), mix4Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array, mix3Array, mix4Array};

                case {'ROCKET_FAC'}
                    % Initialization
                    mix3Array = mix1Array;
                    mix4Array = mix1Array;
                    mix5Array = mix1Array;
                    
                    % Calculations
                    [mix1Array(n), mix2Array(n), mix3Array(n), mix4Array(n), mix5Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i), mix3Array(i), mix4Array(i), mix5Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1), mix3Array(i + 1), mix4Array(i + 1), mix5Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array, mix3Array, mix4Array, mix5Array};

            end

            % Timer
            obj.time = toc(obj.time);

            % Print elapsed time
            printTime(obj);

            % Postprocess all the results with predefined plots
            if obj.FLAG_REPORT
                report(obj, varargout{:});
            end
        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            
            if ~obj.FLAG_TIME
                return
            end

            fprintf('\nElapsed time is %.5f seconds\n', obj.time);
        end

        function ax2 = plot(obj, mixArray1, mixArray2, varargin)
            % Plot results
            %
            % Args:
            %     obj (RocketSolver): RocketSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (initial state)
            %     mixArray2 (Mixture): Array of Mixture objects
            %
            % Optional Args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Examples:
            %     * plot(RocketSolver(), mixArray1, mixArray2, mixArray3);
            %     * plot(RocketSolver(), mixArray1, mixArray2, mixArray3, mixArray4);
            %     * plot(RocketSolver(), mixArray1, mixArray2, mixArray3, mixArray4, mixArray5);
            
            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            additionalMixtures = nargin - 3;
            numPlotProperties = obj.plotConfig.numPlotProperties;
            % indicesToRemoveForInjector = ismember(obj.plotConfig.plotProperties, {'I_sp', 'I_vac'});
            % plotPropertiesInjector = obj.plotConfig.plotProperties(~indicesToRemoveForInjector);
            % plotPropertiesBasisInjector = obj.plotConfig.plotPropertiesBasis(~indicesToRemoveForInjector);
            % numPropertiesInjector = length(plotPropertiesInjector);

            % Check if is a scalar value
            if isscalar(mixArray1)
                ax2 = [];
                return
            end
            
            % Get labels
            switch upper(obj.problemType)
                case {'ROCKET_IAC'}
                    labels = {'Chamber', 'Throat', 'Exit'};
                case {'ROCKET_FAC'}
                    labels = {'Injector', 'Chamber', 'Throat', 'Exit'};

                otherwise
                    if additionalMixtures
                        labels = arrayfun(@(x) sprintf('Mixture %d', x), 1:(additionalMixtures + 1), 'UniformOutput', false);
                    else
                        labels = {''};
                    end
            end
            
            % Plot molar fractions - mixArray2
            ax1 = plotComposition(mixArray2(1), mixArray1, mixArray1(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay, 'y_var', mixArray2, 'title', labels{1});
        
            % Plot properties - mixArray2
            ax2 = plotProperties(repmat({mixArray1(1).rangeName}, 1, numPlotProperties), mixArray1, obj.plotConfig.plotProperties, mixArray2, 'basis', obj.plotConfig.plotPropertiesBasis, 'config', obj.plotConfig);

            % Check if there are additional mixtures
            if ~additionalMixtures
                return
            end

            for i = 1:additionalMixtures
                % Unpack input
                mixArray2 = varargin{i};

                % Plot molar fractions - mixArray_i
                ax1 = plotComposition(mixArray2(1), mixArray1, mixArray1(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay, 'y_var', mixArray2, 'title', labels{i + 1});
            
                % Plot properties - mixArray_i
                ax2 = plotProperties(repmat({mixArray1(1).rangeName}, 1, numPlotProperties), mixArray1, obj.plotConfig.plotProperties, mixArray2, 'basis', obj.plotConfig.plotPropertiesBasis, 'config', obj.plotConfig, 'ax', ax2);
            end

            % Set legends
            legend(ax2.Children(end), labels(1:1+i), 'Interpreter', 'latex', 'FontSize', ax2.Children(end).FontSize);
        end

        function report(obj, mixArray1, mixArray2, varargin)
            % Postprocess all the results with predefined plots
            %
            % Args:
            %     obj (DetonationSolver): DetonationSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (pre-shock state)
            %     mixArray2 (Mixture): Array of Mixture objects (post-shock state)
            %
            % Optional args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Examples:
            %     * report(DetonationSolver(), mixArray1, mixArray2);
            %     * report(DetonationSolver(), mixArray1, mixArray2, mixArray3);

            if nargin > 3
                obj.plot(mixArray1, mixArray2, varargin{:});
            else
                obj.plot(mixArray1, mixArray2);
            end

        end

    end

    methods (Access = private)
        mix = rocketChamberIAC(obj, mix, mix_guess)
        mix3 = rocketThroatIAC(obj, mix2, mix3)
        mix4 = rocketExit(obj, mix2, mix3, mix4, areaRatio, varargin)
        [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1, varargin)
        [mix1, mix2_inj, mix2_c, mix3, mix4] = rocketFAC(obj, mix1, varargin)
    end
    
    methods (Static)
        pressure = rocketGuessThroatIAC(mix)
        log_Pe = rocketGuessExitIAC(mix2, mix3, areaRatio, FLAG_SUBSONIC)
        pressure_inf = rocketGuessInjectorFAC(pressure_inj, areaRatioChamber)
        [mix3, varargout] = rocketParameters(mix2, mix3, varargin)
    end

end
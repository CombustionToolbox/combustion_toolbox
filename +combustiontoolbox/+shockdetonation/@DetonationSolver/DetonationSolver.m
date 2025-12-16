classdef DetonationSolver < handle
    % The :mat:func:`DetonationSolver` class is used to solve detonation waves problems under different
    % flow configurations. The class provides methods to compute the post-detonation state of the gas
    % mixture, including incident and reflected detonations in over-driven and under-driven conditions,
    % oblique detonations, and detonation polar diagrams.
    %
    % The :mat:func:`DetonationSolver` object can be initialized as follows: ::
    %
    %       solver = DetonationSolver('problemType', problemType, ...)
    %
    % Here ``problemType`` represents the acronym of the problem to be solved (see below).
    % Additional optional parameters can be provided to customize the solver's behavior.
    % 
    % Problem types:
    %     * ``DET``: Chapman-Jouguet detonation
    %     * ``DET_R``: Chapman-Jouguet reflected detonation
    %     * ``DET_OVERDRIVEN``: Over-driven detonation
    %     * ``DET_OVERDRIVEN_R``: Reflected over-driven detonation
    %     * ``DET_UNDERDRIVEN``: Under-driven detonation
    %     * ``DET_UNDERDRIVEN_R``: Reflected under-driven detonation
    %     * ``DET_OBLIQUE``: Oblique detonation
    %     * ``DET_OBLIQUE_R``: Oblique reflected detonation 
    %     * ``DET_POLAR``: Detonation polar diagrams
    %     * ``DET_POLAR_R``: Detonation polar diagrams for incident and reflected states
    %
    % See also: :mat:func:`Mixture`, :mat:func:`EquilibriumSolver`, :mat:func:`ShockSolver`, :mat:func:`solve`, :mat:func:`solveArray`, :mat:func:`report`

    properties
        problemType             % Problem type
        equilibriumSolver       % EquilibriumSolver object
        shockSolver             % ShockSolver object
        tol0 = 1e-5             % Tolerance of shocks/detonations kernel
        itMax = 50              % Max number of iterations - shocks and detonations
        machThermo = 2          % Pre-shock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2
        tolOblique = 1e-3       % Tolerance oblique shocks algorithm
        itOblique = 20          % Max number of iterations - oblique shocks
        numPointsPolar = 100    % Number of points to compute shock/detonation polar curves
        tolLimitRR = 1e-4       % Tolerance to calculate the limit of regular reflections
        itLimitRR = 10          % Max number of iterations - limit of regular reflections
        itGuess = 5             % Max number of iterations - guess detonation
        FLAG_RESULTS = true     % Flag to print results
        FLAG_TIME = true        % Flag to print elapsed time
        FLAG_REPORT = false     % Flag to print predefined plots
        FLAG_CACHE = true       % Flag to clear cache after calculations
        time                    % Elapsed time [s]
        plotConfig              % PlotConfig object
    end

    methods
        
        [P, T, STOP] = detonationGuess(obj, mix1)
        [mix1, mix2] = detonationCJ(obj, mix1, varargin)
        [mix1, mix2] = detonationOverdriven(obj, mix1, driveFactor, varargin)
        [mix1, mix2] = detonationUnderdriven(obj, mix1, driveFactor, varargin)
        [mix1, mix2, mix5] = detonationReflected(obj, mix1, mix2, varargin)
        [mix1, mix2_1, mix2_2] = detonationObliqueBeta(obj, mix1, driveFactor, varargin)
        [mix1, mix2_1, mix2_2] = detonationObliqueTheta(obj, mix1, driveFactor, theta, varargin)
        [mix1, mix2] = detonationPolar(obj, mix1, u1)
        [mix1, mix2, mix2_1, mix3] = detonationPolarLimitRR(obj, mix1, u1)

        function obj = DetonationSolver(varargin)
            % Constructor
            defaultProblemType = 'DET';
            defaultCaloricGasModel = combustiontoolbox.core.CaloricGasModel.imperfect;
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('plotConfig', defaultPlotConfig);
            defaultShockSolver = combustiontoolbox.shockdetonation.ShockSolver('plotConfig', defaultPlotConfig);
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;
            
            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'DET', 'DET_R', 'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R', 'DET_OBLIQUE', 'DET_OBLIQUE_R', 'DET_POLAR', 'DET_POLAR_R'})));
            addParameter(p, 'caloricGasModel', defaultCaloricGasModel, @(x) isa(x, 'combustiontoolbox.core.CaloricGasModel'));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addParameter(p, 'shockSolver', defaultShockSolver);
            addParameter(p, 'tol0', obj.tol0, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itMax', obj.itMax, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'machThermo', obj.machThermo, @(x) isnumeric(x) && x >= 1);
            addParameter(p, 'tolOblique', obj.tolOblique, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itOblique', obj.itOblique, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'numPointsPolar', obj.numPointsPolar, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'tolLimitRR', obj.tolLimitRR, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'itLimitRR', obj.itLimitRR, @(x) isnumeric(x) && x > 0);
            addParameter(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_TCHEM_FROZEN', defaultFLAG_TCHEM_FROZEN, @(x) islogical(x))
            addParameter(p, 'FLAG_FROZEN', defaultFLAG_FROZEN, @(x) islogical(x))
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'FLAG_CACHE', obj.FLAG_CACHE, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            addParameter(p, 'tolMoles', defaultEquilibriumSolver.tolMoles, @(x) isnumeric(x) && x > 0);
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.shockSolver = p.Results.shockSolver;
            obj.tol0 = p.Results.tol0;
            obj.itMax = p.Results.itMax;
            obj.machThermo = p.Results.machThermo;
            obj.tolOblique = p.Results.tolOblique;
            obj.itOblique = p.Results.itOblique;
            obj.numPointsPolar = p.Results.numPointsPolar;
            obj.tolLimitRR = p.Results.tolLimitRR;
            obj.itLimitRR = p.Results.itLimitRR;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.FLAG_CACHE = p.Results.FLAG_CACHE;
            obj.plotConfig = p.Results.plotConfig;

            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.caloricGasModel = p.Results.caloricGasModel;
                obj.equilibriumSolver.tolMoles = p.Results.tolMoles;
            end
            
            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
            obj.equilibriumSolver.FLAG_CACHE = false;
            obj.shockSolver.FLAG_RESULTS = false;
            obj.shockSolver.FLAG_TIME = false;
            obj.shockSolver.FLAG_CACHE = false;
            obj.plotConfig.plotProperties{end + 1} = 'uShock';
            obj.plotConfig.plotPropertiesBasis{end + 1} = [];

            % Display warning if deprecated flags are used
            if ~ismember('FLAG_TCHEM_FROZEN', p.UsingDefaults) || ~ismember('FLAG_FROZEN', p.UsingDefaults)
                warning(['The flags ''FLAG_TCHEM_FROZEN'' and ''FLAG_FROZEN'' are deprecated. ', ...
                         'Please use the ''caloricGasModel'' parameter with values from the CaloricGasModel enumeration instead.']);
            
                obj.equilibriumSolver.caloricGasModel = obj.equilibriumSolver.caloricGasModel.fromFlag(p.Results.FLAG_TCHEM_FROZEN, p.Results.FLAG_FROZEN);
            end

        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the DetonationSolver object
            %
            % Args:
            %     obj (DetonationSolver): DetonationSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (DetonationSolver): DetonationSolver object with updated properties
            %
            % Examples:
            %     * set(DetonationSolver(), 'tol0', 1e-6);
            %     * set(DetonationSolver(), 'problemType', 'DET');
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function varargout = solve(obj, mix1, varargin)
            % Solve detonations problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): initial Mixture object
            %
            % Optional Args:
            %     * mixGuess_i (Mixture): Mixture object from previous calculation
            %
            % Returns:
            %     varargout (Mixture): Updated Mixture objects depending on the problem type
            %
            % Examples:
            %     * [mix1, mix2] = solve(DetonationSolver(), mix1);
            %     * [mix1, mix2] = solve(DetonationSolver(), mix1, mix2Guess);
            %     * [mix1, mix2, mix3] = solve(DetonationSolver(), mix1, mix2Guess, mix3Guess);
            
            % Definitions
            driveFactor = mix1.driveFactor;
            beta = mix1.beta;
            theta = mix1.theta;

            switch upper(obj.problemType)
                case 'DET'
                    if nargin > 2
                        [mix1, mix2] = obj.detonationCJ(mix1, varargin{1});
                    else
                        [mix1, mix2] = obj.detonationCJ(mix1);
                    end
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2);
                    end

                    % Set output
                    varargout = {mix1, mix2};

                case 'DET_R'
                    if nargin > 2
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationCJ(mix1, varargin{1});
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2, varargin{2});
                    else
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationCJ(mix1);
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2);
                    end

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix3);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix3};
                
                case 'DET_OVERDRIVEN'
                    if nargin > 2
                        [mix1, mix2] = obj.detonationOverdriven(mix1, driveFactor, varargin{1});
                    else
                        [mix1, mix2] = obj.detonationOverdriven(mix1, driveFactor);
                    end
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
    
                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2);
                    end
    
                    % Set output
                    varargout = {mix1, mix2};
                
                case 'DET_OVERDRIVEN_R'
                    if nargin > 2
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationOverdriven(mix1, driveFactor, varargin{1});
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2, varargin{2});
                    else
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationOverdriven(mix1, driveFactor);
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2);
                    end

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix3);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix3};

                case 'DET_UNDERDRIVEN'
                    if nargin > 2
                        [mix1, mix2] = obj.detonationUnderdriven(mix1, driveFactor, varargin{1});
                    else
                        [mix1, mix2] = obj.detonationUnderdriven(mix1, driveFactor);
                    end
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
    
                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2);
                    end
    
                    % Set output
                    varargout = {mix1, mix2};
                
                case 'DET_UNDERDRIVEN_R'
                    if nargin > 2
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationUnderdriven(mix1, driveFactor, varargin{1});
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2, varargin{2});
                    else
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationUnderdriven(mix1, driveFactor);
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2);
                    end

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix3);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix3};

                case 'DET_OBLIQUE'
                    
                    if isempty(mix1.theta)

                        [mix1, mix2, mix3] = obj.detonationObliqueBeta(mix1, driveFactor, beta);

                        % Set problemType
                        mix1.problemType = obj.problemType;
                        mix2.problemType = obj.problemType;
                        mix3.problemType = obj.problemType;

                        % Print results
                        if obj.FLAG_RESULTS
                            print(mix1, mix2, mix3);
                        end
    
                        % Set output
                        varargout = {mix1, mix2, mix3};
                    else
                        
                        [mix1, mix2, mix3] = obj.detonationObliqueTheta(mix1, driveFactor, theta);
                        
                        % Set problemType
                        mix1.problemType = obj.problemType;
                        mix2.problemType = obj.problemType;
                        mix3.problemType = obj.problemType;
    
                        % Print results
                        if obj.FLAG_RESULTS
                            print(mix1, mix2, mix3);
                        end
    
                        % Set output
                        varargout = {mix1, mix2, mix3};
                    end
                
                case 'DET_POLAR'
                    % Solve problem
                    [mix1, mix2] = obj.detonationPolar(mix1, driveFactor);
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2);
                    end

                    % Set output
                    varargout = {mix1, mix2};

                case 'DET_POLAR_R'
                    % Solve problem
                    [mix1, mix2] = obj.detonationPolar(mix1, driveFactor);
                    [~, mix2_1] = obj.detonationObliqueTheta(mix1, driveFactor, theta);
                    [~, mix3] = obj.detonationPolar(mix2_1, mix2_1.driveFactor);
                    [~, mix3_1, mix3_2] = obj.detonationObliqueTheta(mix2_1, mix2_1.driveFactor, theta);
                
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix2_1, mix3, mix3_1, mix3_2);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix2_1, mix3, mix3_1, mix3_2};

                case 'DET_POLAR_LIMITRR'
                    % Solve problem
                    [mix1, mix2, mix2_1, mix3] = obj.detonationPolarLimitRR(mix1, u1);
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix2_1.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;
                    
                    % Assing values
                    obj.setMixtureShockPolar(mix2_1, mix2);

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix2_1, mix3);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix2_1, mix3};
                    
            end

        end

        function varargout = solveArray(obj, mixArray1, varargin)
            % Solve a set of detonation problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mixArray1 (Mixture): Array of initial Mixture objects
            %
            % Returns:
            %     varargout (Mixture): Updated arrays of Mixture objects depending on the shock problem type
            %
            % Examples:
            %     * [mixArray1, mixArray2] = solveArray(DetonationSolver(), mixArray1);
            %     * [mixArray1, mixArray2, mixArray3] = solveArray(DetonationSolver(), mixArray1);
            
            % Definitions
            n = length(mixArray1);
            problem = obj.problemType;
            
            % Timer
            obj.time = tic;

            % Initialization
            mixArray2 = mixArray1;
            
            % Check conditions
            FLAG_BETA = ~isempty(mixArray1(1).beta);
            FLAG_THETA = ~isempty(mixArray1(1).theta);

            if FLAG_BETA & ~FLAG_THETA
                problem = [problem, '_BETA'];
            elseif ~FLAG_BETA & FLAG_THETA
                problem = [problem, '_THETA'];
            end

            % Calculations
            switch upper(problem)
                case {'DET', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', 'DET_POLAR'}
                    [mixArray1(n), mixArray2(n)] = obj.solve(mixArray1(n));
                    
                    for i = n-1:-1:1
                        [mixArray1(i), mixArray2(i)] = obj.solve(mixArray1(i), mixArray2(i + 1));
                    end

                    % Set output
                    varargout = {mixArray1, mixArray2};

                case {'DET_R', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R', 'DET_OBLIQUE_BETA', 'DET_OBLIQUE_THETA'}
                    % Initialization
                    mixArray3 = mixArray1;
                    
                    % Calculations
                    [mixArray1(n), mixArray2(n), mixArray3(n)] = obj.solve(mixArray1(n));
                    
                    for i = n-1:-1:1
                        [mixArray1(i), mixArray2(i), mixArray3(i)] = obj.solve(mixArray1(i), mixArray2(i + 1), mixArray3(i + 1));
                    end

                    % Set output
                    varargout = {mixArray1, mixArray2, mixArray3};

                case {'DET_POLAR_R_BETA'}
                    % Initialization
                    mixArray3 = mixArray1;
                    mixArray4 = mixArray1;
                    mixArray5 = mixArray1;
                    
                    % Calculations
                    [mixArray1(n), mixArray2(n), mixArray3(n), mixArray4(n), mixArray5(n)] = obj.solve(mixArray1(n));
                    
                    for i = n-1:-1:1
                        [mixArray1(i), mixArray2(i), mixArray3(i), mixArray4(i), mixArray5(i)] = obj.solve(mixArray1(i), mixArray2(i + 1), mixArray3(i + 1), mixArray4(i + 1), mixArray5(i + 1));
                    end

                    % Set output
                    varargout = {mixArray1, mixArray2, mixArray3, mixArray4, mixArray5};

                case {'DET_POLAR_R_THETA'}
                    % Initialization
                    mixArray3 = mixArray1;
                    mixArray4 = mixArray1;
                    mixArray5 = mixArray1;
                    mixArray6 = mixArray1;

                    % Calculations
                    [mixArray1(n), mixArray2(n), mixArray3(n), mixArray4(n), mixArray5(n), mixArray6(n)] = obj.solve(mixArray1(n));
                    
                    for i = n-1:-1:1
                        [mixArray1(i), mixArray2(i), mixArray3(i), mixArray4(i), mixArray5(i), mixArray6(i)] = obj.solve(mixArray1(i), mixArray2(i + 1), mixArray3(i + 1), mixArray4(i + 1), mixArray5(i + 1), mixArray6(i + 1));
                    end

                    % Set output
                    varargout = {mixArray1, mixArray2, mixArray3, mixArray4, mixArray5, mixArray6};

                case {'DET_POLAR_LIMITRR'}
                    % Initialization
                    mixArray3 = mixArray1;
                    mixArray4 = mixArray1;
                    
                    % Calculations
                    [mixArray1(n), mixArray2(n), mixArray3(n), mixArray4(n)] = obj.solve(mixArray1(n));
                    
                    for i = n-1:-1:1
                        [mixArray1(i), mixArray2(i), mixArray3(i), mixArray4(i)] = obj.solve(mixArray1(i), mixArray2(i + 1), mixArray3(i + 1), mixArray4(i + 1));
                    end

                    % Set output
                    varargout = {mixArray1, mixArray2, mixArray3, mixArray4};
            end

            % Timer
            obj.time = toc(obj.time);

            % Print elapsed time
            printTime(obj);

            % Postprocess all the results with predefined plots
            if obj.FLAG_REPORT
                report(obj, varargout{:});
            end


            % Clear cache
            if obj.FLAG_CACHE
                combustiontoolbox.utils.clearCache();
            end

        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            
            if ~obj.FLAG_TIME
                return
            end

            fprintf('\nElapsed time is %.5f seconds\n', obj.time);
        end

        function ax2 = plot(obj, mixArray1, mixArray2, varargin)
            % Plot results
            %
            % Args:
            %     obj (DetonationSolver): DetonationSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (pre-shock state)
            %     mixArray2 (Mixture): Array of Mixture objects (post-shock state)
            %
            % Optional Args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Examples:
            %     * plot(DetonationSolver(), mixArray1, mixArray2);
            %     * plot(DetonationSolver(), mixArray1, mixArray2, mixArray3);
            
            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            additionalMixtures = nargin - 3;
            numPlotProperties = obj.plotConfig.numPlotProperties;

            % Check if is a polar problem
            if contains(obj.problemType, 'POLAR')
                % Plot polars - incident
                plotPolar(mixArray1, mixArray2);
                
                % Check for additional mixtures
                if ~additionalMixtures
                    return
                end

                mixArray2_1 = varargin{1};
                mixArray3 = varargin{2};
                
                plotPolar(mixArray2_1, mixArray3, mixArray2_1, mixArray1);
                return
            end

            % Check if is a scalar value
            if isscalar(mixArray1)
                ax2 = [];
                return
            end
            
            % Get labels
            switch upper(obj.problemType)
                case {'DET_R', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R'}
                    labels = {'Incident', 'Reflected'};
                case {'DET_OBLIQUE'}
                    if isempty(mixArray1(1).theta)
                        labels = {'Underdriven', 'Overdriven'};
                    else
                        labels = {'Weak detonation', 'Strong detonation'};
                    end

                otherwise
                    if additionalMixtures
                        labels = arrayfun(@(x) sprintf('Mixture %d', x), 1:(additionalMixtures + 1), 'UniformOutput', false);
                    else
                        labels = {''};
                    end
            end
            
            % Plot molar fractions - mixArray2
            ax1 = plotComposition(mixArray2(1), mixArray1, mixArray1(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay, 'displaySpecies', obj.plotConfig.displaySpecies, 'y_var', mixArray2, 'title', labels{1});
        
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
                ax1 = plotComposition(mixArray2(1), mixArray1, mixArray1(1).rangeName, 'Xi', 'mintol', obj.plotConfig.mintolDisplay, 'displaySpecies', obj.plotConfig.displaySpecies, 'y_var', mixArray2, 'title', labels{i + 1});
            
                % Plot properties - mixArray_i
                ax2 = plotProperties(repmat({mixArray1(1).rangeName}, 1, numPlotProperties), mixArray1, obj.plotConfig.plotProperties, mixArray2, 'basis', obj.plotConfig.plotPropertiesBasis, 'config', obj.plotConfig, 'ax', ax2);
            end

            % Set legends
            legend(ax2.Children(end), labels, 'Interpreter', 'latex', 'FontSize', ax2.Children(end).FontSize);
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

    methods (Access = private, Static)

        function mix = setMixtureShockPolar(mix, mixPolar)
            % Assign values from the polar curves
            mix.thetaMax = mixPolar.thetaMax;
            mix.betaMax = mixPolar.betaMax;
            mix.thetaSonic = mixPolar.thetaSonic;
            mix.betaSonic = mixPolar.betaSonic;
        end

    end

end
classdef ShockSolver < handle
    % The ShockSolver class is used to solve shock waves problems
    %
    % Problem types:
    %     * SHOCK_I: Incident shock
    %     * SHOCK_R: Reflected shock
    %     * SHOCK_OBLIQUE: Oblique shock
    %     * SHOCK_OBLIQUE_R: Oblique reflected shock 
    %     * SHOCK_POLAR: Shock polar diagrams
    %     * SHOCK_POLAR_R: Shock polar diagrams for incident and reflected states
    %     * SHOCK_POLAR_LIMITRR: Shock polar within the limit of regular reflection
    %
    % Attributes:
    %     problemType (char): Problem type [SHOCK_I, SHOCK_R, SHOCK_OBLIQUE, SHOCK_OBLIQUE_R, SHOCK_POLAR, SHOCK_POLAR_R, SHOCK_POLAR_LIMITRR]
    %     equilibriumSolver (EquilibriumSolver): EquilibriumSolver object
    %     tol0 (float): Tolerance of shocks/detonations kernel
    %     itMax (float): Max number of iterations - shocks and detonations
    %     machThermo (float): Pre-shock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2
    %     tolOblique (float): Tolerance oblique shocks algorithm
    %     itOblique (float): Max number of iterations - oblique shocks
    %     numPointsPolar (float): Number of points to compute shock/detonation polar curves
    %     tolLimitRR (float): Tolerance to calculate the limit of regular reflections
    %     itLimitRR (float): Max number of iterations - limit of regular reflections
    %     FLAG_RESULTS (bool): Flag to print results
    %     FLAG_TIME (bool): Flag to print elapsed time
    %     FLAG_REPORT (bool): Flag to print predefined plots
    %     time (float): Elapsed time
    %
    % Methods:
    %     shockIncident: Solve incident shock
    %     shockReflected: Solve reflected shock
    %     shockObliqueBeta: Solve oblique shock with beta angle
    %     shockObliqueTheta: Solve oblique shock with theta angle
    %     shockObliqueReflectedTheta: Solve oblique reflected shock with theta angle
    %     shockPolar: Solve shock polar diagrams
    %     shockPolarLimitRR: Solve shock polar within the limit of regular reflection
    %     shockIncidentIdeal: Solve incident shock for ideal gas
    %     solve: Solve shock waves problems
    %     solveArray: Solve a set of shock waves problems
    %     printTime: Print execution time
    %     plot: Plot results
    %     report: Postprocess all the results with predefined plots
    %
    % Examples:
    %     * solver = ShockSolver();
    %     * solver = ShockSolver('problemType', 'SHOCK_I');

    properties
        problemType             % Problem type
        equilibriumSolver       % EquilibriumSolver object
        % * Shocks and detonations (CT-SD module)
        tol0 = 1e-5             % Tolerance of shocks/detonations kernel
        itMax = 50              % Max number of iterations - shocks and detonations
        machThermo = 2          % Pre-shock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2
        tolOblique = 1e-3       % Tolerance oblique shocks algorithm
        itOblique = 20          % Max number of iterations - oblique shocks
        numPointsPolar = 100    % Number of points to compute shock/detonation polar curves
        tolLimitRR = 1e-4       % Tolerance to calculate the limit of regular reflections
        itLimitRR = 10          % Max number of iterations - limit of regular reflections
        % * Flags
        FLAG_RESULTS = true     % Flag to print results
        FLAG_TIME = true        % Flag to print elapsed time
        FLAG_REPORT = false     % Flag to print predefined plots
        % * Miscellaneous
        time                    % Elapsed time
        plotConfig              % PlotConfig object
    end

    methods
        
        [mix1, mix2] = shockIncident(obj, mix1, varargin)
        [mix1, mix2, mix5] = shockReflected(obj, mix1, mix2, varargin)
        [mix1, mix2] = shockObliqueBeta(obj, mix1, varargin)
        [mix1, mix2_1, mix2_2] = shockObliqueTheta(obj, mix1, u1, theta, varargin)
        [mix1, mix2, mix5_1, mix5_2] = shockObliqueReflectedTheta(obj, mix1, u2, theta, mix2, varargin)
        [mix1, mix2] = shockPolar(obj, mix1, u1)
        [mix1, mix2, mix2_1, mix3] = shockPolarLimitRR(obj, mix1, u1)
        [R, P, T, Gammas, beta, M1] = shockIncidentIdeal(obj, gamma, M1)

        function obj = ShockSolver(varargin)
            % Constructor
            defaultProblemType = 'SHOCK_I';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;
            
            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'SHOCK_I', 'SHOCK_R', 'SHOCK_OBLIQUE', 'SHOCK_OBLIQUE_R', 'SHOCK_POLAR', 'SHOCK_POLAR_R', 'SHOCK_POLAR_LIMITRR'})));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
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
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
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
            obj.plotConfig = p.Results.plotConfig;

            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
                obj.equilibriumSolver.FLAG_FROZEN = p.Results.FLAG_FROZEN;
            end
            
            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the ShockSolver object
            %
            % Args:
            %     obj (ShockSolver): ShockSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (ShockSolver): ShockSolver object with updated properties
            %
            % Examples:
            %     * set(ShockSolver(), 'tol0', 1e-6);
            %     * set(ShockSolver(), 'problemType', 'SHOCK_I');
            
            varargin = [property, value, varargin];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function varargout = solve(obj, mix1, varargin)
            % Solve shock waves problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1 (Mixture): initial Mixture object
            %
            % Optional Args:
            %     * mixGuess_i (Mixture): Mixture object from previous calculation
            %
            % Returns:
            %     varargout: Updated Mixture objects depending on the problem type
            %
            % Examples:
            %     * [mix1, mix2] = solve(ShockSolver(), mix1);
            %     * [mix1, mix2] = solve(ShockSolver(), mix1, mix2Guess);
            %     * [mix1, mix2, mix3] = solve(ShockSolver(), mix1, mix2Guess, mix3Guess);
            
            % Definitions
            u1 = mix1.u;
            beta = mix1.beta;
            theta = mix1.theta;

            switch upper(obj.problemType)
                case 'SHOCK_I'
                    if nargin > 2
                        [mix1, mix2] = obj.shockIncident(mix1, u1, varargin{1});
                    else
                        [mix1, mix2] = obj.shockIncident(mix1, u1);
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

                case 'SHOCK_R'
                    if nargin > 2
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.shockIncident(mix1, u1, varargin{1});
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.shockReflected(mix1, mix2, varargin{2});
                    else
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.shockIncident(mix1, u1);
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.shockReflected(mix1, mix2);
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

                case 'SHOCK_OBLIQUE'
                    
                    if isempty(mix1.theta)
                        
                        % The initialization using the previous solution is prone to convergence issues.
                        % if nargin > 3
                        %     [mix1, mix2] = obj.shockObliqueBeta(mix1, u1, beta, varargin{1});
                        % else
                        %     [mix1, mix2] = obj.shockObliqueBeta(mix1, u1, beta);
                        % end

                        [mix1, mix2] = obj.shockObliqueBeta(mix1, u1, beta);

                        % Set problemType
                        mix1.problemType = obj.problemType;
                        mix2.problemType = obj.problemType;
    
                        % Print results
                        if obj.FLAG_RESULTS
                            print(mix1, mix2);
                        end
    
                        % Set output
                        varargout = {mix1, mix2};
                    else
                        
                        % The initialization using the previous solution is prone to convergence issues.
                        % if nargin > 3
                        %     [mix1, mix2, mix3] = obj.shockObliqueTheta(mix1, u1, theta, varargin{:});
                        % else
                        %     [mix1, mix2, mix3] = obj.shockObliqueTheta(mix1, u1, theta);
                        % end
                        
                        [mix1, mix2, mix3] = obj.shockObliqueTheta(mix1, u1, theta);
                        
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

                case 'SHOCK_OBLIQUE_R'
                    
                    if isempty(mix1.theta)
                        % Compute incident shock
                        [mix1, mix2] = obj.shockObliqueBeta(mix1, u1, beta);
                    else
                        % Compute incident shock
                        [mix1, mix2, ~] = obj.shockObliqueTheta(mix1, u1, theta);
                    end

                    % Compute reflected shock
                    [mix1, mix2, mix3_1, mix3_2] = obj.shockObliqueReflectedTheta(mix1, mix2.u, mix2.theta, mix2);

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix3_1.problemType = obj.problemType;
                    mix3_2.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2, mix3_1, mix3_2);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix3_1, mix3_2};
                
                case 'SHOCK_POLAR'
                    % Solve problem
                    [mix1, mix2] = obj.shockPolar(mix1, u1);
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2);
                    end

                    % Set output
                    varargout = {mix1, mix2};

                case 'SHOCK_POLAR_R'
                    % Solve problem
                    [mix1, mix2] = obj.shockPolar(mix1, u1);
                    [~, mix2_1] = obj.shockObliqueTheta(mix1, u1, theta);
                    [~, mix3] = obj.shockPolar(mix2_1, mix2_1.u);
                    [~, mix3_1, mix3_2] = obj.shockObliqueTheta(mix2_1, mix2_1.u, theta);
                
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    
                    % Assing values
                    obj.setMixtureShockPolar(mix2_1, mix2);
                    obj.setMixtureShockPolar(mix3_1, mix3);
                    obj.setMixtureShockPolar(mix3_2, mix3);

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2_1, mix3_1, mix3_2);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix2_1, mix3, mix3_1, mix3_2};

                case 'SHOCK_POLAR_LIMITRR'
                    % Solve problem
                    [mix1, mix2, mix2_1, mix3] = obj.shockPolarLimitRR(mix1, u1);
                    
                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2.problemType = obj.problemType;
                    mix2_1.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;
                    
                    % Assing values
                    obj.setMixtureShockPolar(mix2_1, mix2);

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2_1, mix3);
                    end

                    % Set output
                    varargout = {mix1, mix2, mix2_1, mix3};
            end

        end

        function varargout = solveArray(obj, mix1Array, varargin)
            % Solve a set of shock waves problems
            %
            % Args:
            %     obj (EquilibriumSolver): EquilibriumSolver object
            %     mix1Array (Mixture): Array of initial Mixture objects
            %
            % Returns:
            %     varargout: Updated arrays of Mixture objects depending on the shock problem type
            %
            % Examples:
            %     * [mix1Array, mix2Array] = solveArray(ShockSolver(), mix1Array);
            %     * [mix1Array, mix2Array, mix3Array] = solveArray(ShockSolver(), mix1Array);

            % Definitions
            n = length(mix1Array);
            problem = obj.problemType;
            
            % Timer
            obj.time = tic;

            % Initialization
            mix2Array = mix1Array;
            
            % Check conditions
            FLAG_BETA = ~isempty(mix1Array(1).beta);
            FLAG_THETA = ~isempty(mix1Array(1).theta);

            if FLAG_BETA & ~FLAG_THETA
                problem = [problem, '_BETA'];
            elseif ~FLAG_BETA & FLAG_THETA
                problem = [problem, '_THETA'];
            end

            % Calculations
            switch upper(problem)
                case {'SHOCK_I', 'SHOCK_OBLIQUE_BETA', 'SHOCK_POLAR'}
                    [mix1Array(n), mix2Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array};

                case {'SHOCK_R', 'SHOCK_OBLIQUE_THETA'}
                    % Initialization
                    mix3Array = mix1Array;
                    
                    % Calculations
                    [mix1Array(n), mix2Array(n), mix3Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i), mix3Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1), mix3Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array, mix3Array};

                case {'SHOCK_POLAR_R_BETA'}
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

                case {'SHOCK_POLAR_R_THETA'}
                    % Initialization
                    mix3Array = mix1Array;
                    mix4Array = mix1Array;
                    mix5Array = mix1Array;
                    mix6Array = mix1Array;

                    % Calculations
                    [mix1Array(n), mix2Array(n), mix3Array(n), mix4Array(n), mix5Array(n), mix6Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i), mix3Array(i), mix4Array(i), mix5Array(i), mix6Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1), mix3Array(i + 1), mix4Array(i + 1), mix5Array(i + 1), mix6Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array, mix3Array, mix4Array, mix5Array, mix6Array};

                case {'SHOCK_OBLIQUE_R_BETA', 'SHOCK_OBLIQUE_R_THETA', 'SHOCK_POLAR_LIMITRR'}
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
            %     obj (ShockSolver): ShockSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (pre-shock state)
            %     mixArray2 (Mixture): Array of Mixture objects (post-shock state)
            %
            % Optional Args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Examples:
            %     * plot(ShockSolver(), mixArray1, mixArray2);
            %     * plot(ShockSolver(), mixArray1, mixArray2, mixArray3);
            
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
                case 'SHOCK_R'
                    labels = {'Incident', 'Reflected'};
                case 'SHOCK_OBLIQUE'
                    if isempty(mixArray1(1).theta)
                        labels = {''};
                    else
                        labels = {'Weak shock', 'Strong shock'};
                    end

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

            % Plot Hugoniot curve
            % ax3 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'color', 'auto');
            
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

                % Plot Hugoniot curve
                % ax3 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'color', 'auto', 'ax', ax3);
            end

            % Set legends
            legend(ax2.Children(end), labels, 'Interpreter', 'latex', 'FontSize', ax2.Children(end).FontSize);
            % legend(ax3, legendName, 'Interpreter', 'latex', 'FontSize', ax3.FontSize);
        end

        function report(obj, mixArray1, mixArray2, varargin)
            % Postprocess all the results with predefined plots
            %
            % Args:
            %     obj (ShockSolver): ShockSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (pre-shock state)
            %     mixArray2 (Mixture): Array of Mixture objects (post-shock state)
            %
            % Optional args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Examples:
            %     * report(ShockSolver(), mixArray1, mixArray2);
            %     * report(ShockSolver(), mixArray1, mixArray2, mixArray3);

            if nargin > 3
                obj.plot(mixArray1, mixArray2, varargin{:});
            else
                obj.plot(mixArray1, mixArray2);
            end

        end

    end
    
    methods (Access = public, Static)

        function Gammas = computeGammas(u2, rho2, p2)
            % Compute slope of the Hugoniot curve
            %
            % Args:
            %     obj (obj): 
            %     u2 (float): Post-shock velocity [m/s]
            %     rho2 (float): Post-shock density [kg/m3]
            %     p2 (float): Post-shock pressure [bar]
            %
            % Returns:
            %     Gammas (float): Slope of the Hugoniot curve [-]
            
            p2 = convert_bar_to_Pa(p2);
            
            Gammas =  u2(1:end-1).^2 .* combustiontoolbox.utils.math.computeFirstDerivative(rho2, p2);
        end

        function Gammas = computeGammasFrozen(M1, R, P)
            % Compute slope of the Hugoniot curve for thermochemically frozen air
            %
            % Args:
            %     M1 (float): Pre-shock Mach number [-]
            %     R (float): Density ratio [-]
            %     P (float): Pressure ratio [-]
            %
            % Returns:
            %     Gammas (float): Slope of the Hugoniot curve [-]
        
            Gammas =  7/5 * (M1(1:end-1).^2 ./ R(1:end-1).^2) .* combustiontoolbox.utils.math.computeFirstDerivative(P, R).^(-1);
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
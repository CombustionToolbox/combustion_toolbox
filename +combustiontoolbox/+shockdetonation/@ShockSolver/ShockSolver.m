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
    %
    % Examples:
    %     * solver = ShockSolver();
    %     * solver = ShockSolver('problemType', 'SHOCK_I');

    properties
        problemType             % Problem type
        equilibriumSolver       % EquilibriumSolver
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
        % * Miscellaneous
        time
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

            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
                obj.equilibriumSolver.FLAG_FROZEN = p.Results.FLAG_FROZEN;
            end
            
            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
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
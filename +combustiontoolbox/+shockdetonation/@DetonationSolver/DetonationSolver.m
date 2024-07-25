classdef DetonationSolver

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
        itGuess = 5             % Max number of iterations - guess detonation
        % * Flags
        FLAG_RESULTS = true     % Flag to print results
        FLAG_TIME = true        % Flag to print elapsed time
        % * Miscellaneous
        time
    end

    methods
        
        [mix1, mix2] = detonationCJ(obj, mix1, varargin)
        [P, T, STOP] = detonationGuess(obj, mix1)
        [mix1, mix2] = detonationOverdriven(obj, mix1, driveFactor, varargin)
        [mix1, mix2] = detonationUnderdriven(obj, mix1, driveFactor, varargin)
        [mix1, mix2, mix5] = detonationReflected(obj, mix1, mix2, varargin)
        [mix1, mix2] = detonationObliqueBeta(obj, mix1, varargin)
        [mix1, mix2_1, mix2_2] = detonationObliqueTheta(obj, mix1, u1, theta, varargin)
        [mix1, mix2] = detonationPolar(obj, mix1, u1)
        [mix1, mix2, mix2_1, mix3] = detonationPolarLimitRR(obj, mix1, u1)

        function obj = DetonationSolver(varargin)
            % Constructor
            defaultProblemType = 'DET';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;
            
            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'DET', 'DET_R', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', 'DET_OBLIQUE', 'DET_OBLIQUE_R', 'DET_POLAR', 'DET_POLAR_R', 'DET_POLAR_LIMITRR'})));
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
            % Solve detonations problems
            
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
                        [mix1, mix2] = obj.detonationCJ(mix1, u1, varargin{1});
                        % Calculate post-shock state (5)
                        [mix1, mix2, mix3] = obj.detonationReflected(mix1, mix2, varargin{2});
                    else
                        % Calculate post-shock state (2)
                        [mix1, mix2] = obj.detonationCJ(mix1, u1);
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

                case 'DET_OBLIQUE'
                    
                    if isempty(mix1.theta)

                        [mix1, mix2] = obj.detonationObliqueBeta(mix1, u1, beta);

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
                        
                        [mix1, mix2, mix3] = obj.detonationObliqueTheta(mix1, u1, theta);
                        
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
                    [mix1, mix2] = obj.detonationPolar(mix1, u1);
                    
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
                    [mix1, mix2] = obj.detonationPolar(mix1, u1);
                    [~, mix2_1] = obj.detonationObliqueTheta(mix1, u1, theta);
                    [~, mix3] = obj.detonationPolar(mix2_1, mix2_1.u);
                    [~, mix3_1, mix3_2] = obj.detonationObliqueTheta(mix2_1, mix2_1.u, theta);
                
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

        function varargout = solveArray(obj, mix1Array, varargin)
            % Solve a set of detonations problems
            
            % Definitions
            n = length(mix1Array);
            problem = obj.problemType;

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
                case {'DET', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', 'DET_OBLIQUE_BETA', 'DET_POLAR'}
                    [mix1Array(n), mix2Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array};

                case {'DET_R', 'DET_OBLIQUE_THETA'}
                    % Initialization
                    mix3Array = mix1Array;
                    
                    % Calculations
                    [mix1Array(n), mix2Array(n), mix3Array(n)] = obj.solve(mix1Array(n));
                    
                    for i = n-1:-1:1
                        [mix1Array(i), mix2Array(i), mix3Array(i)] = obj.solve(mix1Array(i), mix2Array(i + 1), mix3Array(i + 1));
                    end

                    % Set output
                    varargout = {mix1Array, mix2Array, mix3Array};

                case {'DET_POLAR_R_BETA'}
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

                case {'DET_POLAR_R_THETA'}
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

                case {'DET_POLAR_LIMITRR'}
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
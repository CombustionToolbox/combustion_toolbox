classdef ShockSolver < handle

    properties
        problemType             % Problem type
        equilibriumSolver       % EquilibriumSolver
        % * Shocks and detonations (CT-SD module)
        tol_shocks = 1e-5       % Tolerance of shocks/detonations kernel
        it_shocks = 50          % Max number of iterations - shocks and detonations
        Mach_thermo = 2         % Pre-shock Mach number above which T2_guess will be computed considering h2 = h1 + u1^2 / 2
        tol_oblique = 1e-3      % Tolerance oblique shocks algorithm
        it_oblique = 20         % Max number of iterations - oblique shocks
        N_points_polar = 100    % Number of points to compute shock/detonation polar curves
        tol_limitRR = 1e-4      % Tolerance to calculate the limit of regular reflections
        it_limitRR = 10         % Max number of iterations - limit of regular reflections
        it_guess_det = 5        % Max number of iterations - guess detonation
        % Miscellaneous
        FLAG_RESULTS = true     % Flag to print results
    end

    methods

        function obj = ShockSolver(varargin)
            % Constructor
            defaultProblemType = 'SHOCK_I';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'SHOCK_I', 'SHOCK_R', 'SHOCK_OBLIQUE', 'SHOCK_POLAR'})));
            addOptional(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addOptional(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;

            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
        end

        function varargout = solve(obj, mix1, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties
            
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
                    
            end

        end

        function varargout = solveArray(obj, mix1Array, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties
            %
            %
            %
            
            % Definitions
            n = length(mix1Array);
            problem = obj.problemType;

            % Initialization
            mix2Array = mix1Array;
            
            % Check conditions
            FLAG_BETA = ~isempty(mix1Array(1).beta);
            FLAG_THETA = ~isempty(mix1Array(1).theta);

            if FLAG_BETA & ~FLAG_THETA
                problem = 'SHOCK_OBLIQUE_BETA';
            elseif ~FLAG_BETA & FLAG_THETA
                problem = 'SHOCK_OBLIQUE_THETA';
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
            end

        end

    end

    methods (Access = private)
        
        [mix1, mix2] = shockIncident(obj, mix1, varargin)
        [mix1, mix2, mix5] = shockReflected(obj, mix1, mix2, varargin)
        [mix1, mix2] = shockObliqueBeta(obj, mix1, varargin)
        [mix1, mix2_1, mix2_2] = shockObliqueTheta(obj, mix1, u1, theta, varargin);
        [mix1, mix2] = shockPolar(obj, mix1, u1)

    end

end
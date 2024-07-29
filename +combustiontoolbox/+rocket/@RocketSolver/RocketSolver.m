classdef RocketSolver

    properties
        problemType           % Problem type
        equilibriumSolver     % EquilibriumSolver
        % * Rocket propellant performance (CT-ROCKET module)
        tol0 = 1e-4;          % Tolerance rocket performance
        itMax = 10;           % Max number of iterations - rocket performance
        % * Flags
        FLAG_SUBSONIC = false % Flag to indicate subsonic Area ratio (CT-ROCKET)
        FLAG_RESULTS = true   % Flag to print results
        FLAG_TIME = true      % Flag to print elapsed time
        % * Miscellaneous
        time
    end

    methods

        function obj = RocketSolver(varargin)
            % Constructor
            defaultProblemType = 'ROCKET_IAC';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'ROCKET_IAC', 'ROCKET_FAC'})));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addParameter(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;

            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
        end

        function varargout = solve(obj, mix1, varargin)
            % Solve rocket problems
            
            switch upper(obj.problemType)
                case 'ROCKET_IAC'
                    % Solve rocket problem with Infinite Area Chamber (IAC) model
                    [mix1, mix2_c, mix3, mix4] = rocketIAC(obj, mix1);

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2_c.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;

                    % Print results
                    if obj.FLAG_RESULTS
                        print(mix1, mix2_c, mix3, mix4);
                    end

                    % Set output
                    varargout = {mix1, mix2_c, mix3, mix4};
                case 'ROCKET_FAC'
                    % Solve rocket problem considering an Finite Area Chamber (FAC) model
                    [mix1, mix2_inj, mix2_c, mix3, mix4] = rocketFAC(obj, mix1);

                    % Set problemType
                    mix1.problemType = obj.problemType;
                    mix2_inj.problemType = obj.problemType;
                    mix2_c.problemType = obj.problemType;
                    mix3.problemType = obj.problemType;

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
            % Solve a set of detonations problems
            
            % Definitions
            n = length(mix1Array);
            problem = obj.problemType;

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

        end

    end

    methods (Access = private)
        mix = computeChamberIAC(obj, mix)
        mix3 = computeThroatIAC(obj, mix2, mix3)
        mix4 = computeExit(obj, mix2, mix3, mix4, areaRatio, varargin)
        [mix3, varargout] = rocketParameters(obj, mix2, mix3, varargin)
    end

end
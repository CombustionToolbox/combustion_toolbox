classdef RocketSolver
    % RocketSolver class to solve rocket problems
    %
    % Problem types:
    %   * ROCKET_IAC: Get the performance of a rocket using the Infinite Area Chamber (IAC) model
    %   * ROCKET_FAC: Get the performance of a rocket using the Finite Area Chamber (FAC) model
    %
    % Attributes:
    %     problemType (char): Problem type [ROCKET_IAC or ROCKET_FAC]
    %     equilibriumSolver (EquilibriumSolver): EquilibriumSolver object
    %     tol0 (float): Tolerance rocket performance
    %     itMax (float): Max number of iterations - rocket performance
    %     FLAG_SUBSONIC (bool): Flag to indicate subsonic Area ratio
    %     FLAG_RESULTS (bool): Flag to print results
    %     FLAG_TIME (bool): Flag to print elapsed time
    %     time (float): Elapsed time
    %
    % Methods:
    %     solve: Solve rocket problems
    %     solveArray: Solve a set of rocket problems
    %     printTime: Print execution time
    %
    % Examples:
    %     * solver = RocketSolver();
    %     * solver = RocketSolver('problemType', 'ROCKET_IAC');
    %     * solver = RocketSolver('problemType', 'ROCKET_FAC');
    %     * solver = RocketSolver('problemType', 'ROCKET_FAC', 'tolMoles', 1e-30);
    %     * solver = RocketSolver('problemType', 'ROCKET_FAC', 'tolMoles', 1e-30, 'FLAG_FROZEN', true);

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
        % * Miscellaneous
        time
    end

    methods

        function obj = RocketSolver(varargin)
            % Constructor
            defaultProblemType = 'ROCKET_IAC';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();
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
            addParameter(p, 'tolMoles', defaultEquilibriumSolver.tolMoles, @(x) isnumeric(x) && x > 0);
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            
            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
                obj.equilibriumSolver.FLAG_FROZEN = p.Results.FLAG_FROZEN;
                obj.equilibriumSolver.tolMoles = p.Results.tolMoles;
            end

            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
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
            %     varargout: Updated Mixture objects depending on the problem type
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
            %     varargout: Updated arrays of Mixture objects depending on the shock problem type
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
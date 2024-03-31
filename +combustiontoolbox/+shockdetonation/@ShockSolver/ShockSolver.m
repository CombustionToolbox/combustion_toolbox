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
    end

    methods

        function obj = ShockSolver(varargin)
            % Constructor
            defaultProblemType = 'SHOCK_I';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();

            % Parse input arguments
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'SHOCK_I'})));
            addOptional(p, 'equilibriumSolver', defaultEquilibriumSolver);
            parse(p, varargin{:});

            % Set properties
            obj.problemType = upper(p.Results.problemType);
            obj.equilibriumSolver = p.Results.equilibriumSolver;
        end

        function [mix1, mix2] = solve(obj, mix1, u1, varargin)
            % Obtain chemical equilibrium composition and thermodynamic properties
            switch obj.problemType
                case 'SHOCK_I'
                    if nargin > 3
                        [mix1, mix2] = obj.shockIncident(mix1, u1, varargin{1});
                        return
                    end
            
                    [mix1, mix2] = obj.shockIncident(mix1, u1);
            end

        end

    end

    methods (Access = private)
        
        [mix1, mix2] = shockIncident(obj, mix1, varargin)

    end

end
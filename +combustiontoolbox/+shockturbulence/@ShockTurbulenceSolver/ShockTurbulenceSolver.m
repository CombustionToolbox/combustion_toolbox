classdef ShockTurbulenceSolver < handle
    % The :mat:func:`ShockTurbulenceSolver` class characterizes turbulence 
    % amplification across a shock wave interacting with weak turbulence using 
    % linear interaction analysis.
    %
    % Uptream turbulence can be comprised of the following type of disturbances:
    % 	* vortical (shockTurbulenceModelVortical)
    %   * acoustic (shockTurbulenceModelAcoustic)
    %   * vortical + entropic (shockTurbulenceModelVorticalEntropic)
    %   * vortical + entropic + acoustic (shockTurbulenceModelCompressible)
    % 
    % These models are based on our previous theoretical work [1]
    % and have been extended to multi-component mixtures [2-4] using the
    % Combustion Toolbox [5].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [3] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [4] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [5] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.

    properties
	    problemType                % Problem type
        equilibriumSolver          % EquilibriumSolver object
        shockSolver                % ShockSolver object
        jumpConditionsSolver       % JumpConditionsSolver object
        shockTurbulenceModel       % ShockTurbulenceModel object
        FLAG_RESULTS = false;      % Flag to show results in the command window
        FLAG_INTERPOLATE = true;   % Flag to interpolate data in a smaller grid
        FLAG_TIME = true           % Flag to print elapsed time
        FLAG_REPORT = false        % Flag to print predefined plots
        time                       % Elapsed time
        plotConfig                 % PlotConfig object
    end

    properties (Constant, Access = private)
        FLAG_PAPER = false;        % Flag to compute Gammas_i as in Cuadra2024b
    end

    methods

        function obj = ShockTurbulenceSolver(varargin)
            % Constructor of the ShockTurbulenceSolver class

            % Default values
            defaultProblemType = 'VORTICAL';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('FLAG_FAST', false); % FLAG_FAST is set to false to reduce numerical error
            defaultShockSolver = combustiontoolbox.shockdetonation.ShockSolver('equilibriumSolver', defaultEquilibriumSolver, 'FLAG_RESULTS', false);
            defaultJumpConditionsSolver = combustiontoolbox.shockdetonation.JumpConditionsSolver('equilibriumSolver', defaultEquilibriumSolver, 'shockSolver', defaultShockSolver,'FLAG_RESULTS', false);
            defaultShockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical();
            defaultCaloricGasModel = combustiontoolbox.core.CaloricGasModel.imperfect;
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultPlotConfig.plotProperties = {'K', 'R11', 'RTT', 'Ka', 'Kr', 'enstrophy'};

            % Parse inputs
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'compressible', 'vortical', 'acoustic', 'vortical_entropic'})));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addParameter(p, 'shockSolver', defaultShockSolver, @(x) isa(x, 'combustiontoolbox.shockdetonation.ShockSolver'));
            addParameter(p, 'jumpConditionsSolver', defaultJumpConditionsSolver, @(x) isa(x, 'combustiontoolbox.shockdetonation.JumpConditionsSolver'));
            addParameter(p, 'shockTurbulenceModel', defaultShockTurbulenceModel, @(x) isa(x, 'combustiontoolbox.shockturbulence.ShockTurbulenceModel'));
            addParameter(p, 'caloricGasModel', defaultCaloricGasModel, @(x) isa(x, 'combustiontoolbox.core.CaloricGasModel'));
            addParameter(p, 'FLAG_INTERPOLATE', obj.FLAG_INTERPOLATE, @islogical);
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'FLAG_TCHEM_FROZEN', defaultFLAG_TCHEM_FROZEN, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_FROZEN', defaultFLAG_FROZEN, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = p.Results.problemType;
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.shockSolver = p.Results.shockSolver;
            obj.jumpConditionsSolver = p.Results.jumpConditionsSolver;
            obj.shockTurbulenceModel = p.Results.shockTurbulenceModel;
            obj.FLAG_INTERPOLATE = p.Results.FLAG_INTERPOLATE;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;

            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.caloricGasModel = p.Results.caloricGasModel;
            end

            % Display warning if deprecated flags are used
            if ~ismember('FLAG_TCHEM_FROZEN', p.UsingDefaults) || ~ismember('FLAG_FROZEN', p.UsingDefaults)
                warning(['The flags ''FLAG_TCHEM_FROZEN'' and ''FLAG_FROZEN'' are deprecated. ', ...
                         'Please use the ''caloricGasModel'' parameter with values from the CaloricGasModel enumeration instead.']);
            
                obj.equilibriumSolver.caloricGasModel = obj.equilibriumSolver.caloricGasModel.fromFlag(p.Results.FLAG_TCHEM_FROZEN, p.Results.FLAG_FROZEN);
            end

            % Assign equilibriumSolver to shockSolver and jumpConditionsSolver
            obj.shockSolver.equilibriumSolver = obj.equilibriumSolver;
            obj.jumpConditionsSolver.equilibriumSolver = obj.equilibriumSolver;

            % Assign shockSolver to jumpConditionsSolver
            obj.jumpConditionsSolver.shockSolver = obj.shockSolver;

            % if jumpConditionsSolver is not default
            if ~sum(contains(p.UsingDefaults, 'jumpConditionsSolver'))
                obj.jumpConditionsSolver.FLAG_PAPER = obj.FLAG_PAPER;
            end

            % If problemType is different from the default, set the corresponding shockTurbulenceModel
            if ~sum(contains(p.UsingDefaults, 'shockTurbulenceModel'))
                return
            end
            
            switch lower(obj.problemType)
                case {'vortical'}
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical();
                case {'vortical_entropic'}
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic();
                case {'acoustic'}
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic();
                case {'compressible'}
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible();
                otherwise
                    error('Unknown problem type: %s', obj.problemType);
            end

        end

        function results = solve(obj, mixArray, varargin)
            % Solve shock waves problems
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     mixArray (Mixture): Initial Mixture objects
            %
            % Returns:
            %     results (struct): Struct with results from LIA
            %
            % Example:
            %     * results = solve(ShockTurbulenceSolver(), mixArray);
            
            % Default
            eta = 0;
            chi = 0;
            jumpConditions = [];

            % Unpack additional inputs
            for i = 1:2:nargin-2

                switch lower(varargin{i})
                    case {'jumpconditions'}
                        jumpConditions = varargin{i + 1};
                    case {'compressibility', 'eta'}
                        eta = varargin{i + 1};
                    case {'vortical_entropic', 'chi'}
                        chi = varargin{i + 1};
                end

            end

            % Timer
 	        obj.time = tic;

            % Compute jump conditions
            if isempty(jumpConditions)
                jumpConditions = obj.jumpConditionsSolver.solve(mixArray, varargin{:});
            end

            % Get jump conditions
            Gammas1 = jumpConditions.Gammas1;
            Gammas2 = jumpConditions.Gammas2;
            Gammas3 = jumpConditions.Gammas3;
            Rratio = jumpConditions.Rratio;
            Pratio = jumpConditions.Pratio;
            Tratio = jumpConditions.Tratio;
            M1 = jumpConditions.M1;
            M2 = jumpConditions.M2;
            beta = jumpConditions.beta;
	    
	        % Solve problem
            switch lower(obj.problemType)
                case {'vortical'}
                    results = obj.shockTurbulenceModel.getAverages(Rratio, M2, Gammas2);
                case {'vortical_entropic'}
                    results = obj.shockTurbulenceModel.getAverages(Rratio, M2, Gammas2, Gammas1, beta, chi);
                case {'acoustic'}
                    results = obj.shockTurbulenceModel.getAverages(Rratio, M2, Gammas2, Gammas1, Gammas3, beta);
                case {'compressible'}
                    results = obj.shockTurbulenceModel.getAverages(Rratio, M2, Gammas2, Gammas1, Gammas3, beta, eta, chi);
            end

            % Add jump conditions
            results.Gammas1 = Gammas1;
            results.Gammas2 = Gammas2;
            results.Gammas3 = Gammas3;
            results.Rratio = Rratio;
            results.Pratio = Pratio;
            results.Tratio = Tratio;
            results.M1 = M1;
            results.M2 = M2;
            results.beta = beta;

            % Compute Kolmogorov length scale ratio across the shock
            results.kolmogorovLengthRatio = obj.getKolmogorovLength(results);

            % Timer
            obj.time = toc(obj.time);

            % Print elapsed time
            printTime(obj);

            % Postprocess all the results with predefined plots
            if obj.FLAG_REPORT
                report(obj, results);
            end

        end

        function printTime(obj)
            % Print execution time
            %
            % Args:
            %     obj (ShockTurbulenceSolver()): Object of the class ShockTurbulenceSolver
            
            if ~obj.FLAG_TIME
                return
            end

            fprintf('\nElapsed time is %.5f seconds\n', obj.time);
        end

        function ax1 = plot(obj, results)
            % Plot results
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     results (struct): Results from the Linear Interaction Analysis
            %
            % Example:
            %     * plot(ShockTurbulenceSolver(), results);
            
            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            numPlotProperties = obj.plotConfig.numPlotProperties;

            % Check if is a scalar value
            if isscalar(results.K)
                ax1 = [];
                return
            end
            
            % Plot properties
            ax1 = plotProperties(repmat({'M1'}, 1, numPlotProperties), results, obj.plotConfig.plotProperties, results, 'config', obj.plotConfig);
        end

        function report(obj, results)
            % Postprocess all the results with predefined plots
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     results (struct): Results from the Linear Interaction Analysis
            %
            % Example:
            %     * report(ShockTurbulenceSolver(), results);

            obj.plot(results);
        end

        function kolmogorovLengthRatio = getKolmogorovLength(obj, results)
            % Estimate Kolmogorov length scale ratio across the shock
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     results (struct): Struct with results from LIA
            %
            % Returns:
            %     kolmogorovLengthRatio (float): Kolmogorov length scale ratio
            %
            % Example:
            %     kolmogorovLengthRatio = getKolmogorovLength(ShockTurbulenceSolver(), results);

            % Definitions
            Rratio = results.Rratio;
            TRatio = results.Tratio;
            enstrophyRatio = results.enstrophy;

            % Kolmogorov length scale ratio 
            kolmogorovLengthRatio = Rratio.^(-1/2) .* TRatio.^(3/8) .* enstrophyRatio.^(-1/4);
        end

    end

end

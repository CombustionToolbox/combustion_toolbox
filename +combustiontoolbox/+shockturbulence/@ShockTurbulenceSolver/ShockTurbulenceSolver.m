classdef ShockTurbulenceSolver < handle
    % The :mat:class:`ShockTurbulenceSolver` class characterizes turbulence 
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
        FLAG_TCHEM_FROZEN = false; % Flag to indicate if the thermodynamic properties are thermochemically frozen (calorically perfect gas)
        FLAG_FROZEN = false;       % Flag to indicate if frozen chemistry (calorically imperfect gas with frozen chemistry)
        FLAG_INTERPOLATE = true;   % Flag to interpolate data in a smaller grid
        FLAG_PAPER = false;        % Flag to compute Gammas_i as in Cuadra2024b
        FLAG_TIME = true           % Flag to print elapsed time
        FLAG_REPORT = false        % Flag to print predefined plots
        time                       % Elapsed time
        plotConfig                 % PlotConfig object
    end

    methods

        function obj = ShockTurbulenceSolver(varargin)
            % Constructor of the ShockTurbulenceSolver class

            % Default values
            defaultProblemType = 'VORTICAL';
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('FLAG_FAST', false, 'FLAG_TCHEM_FROZEN', obj.FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', obj.FLAG_FROZEN); % FLAG_FAST is set to false to reduce numerical error
            defaultShockSolver = combustiontoolbox.shockdetonation.ShockSolver('equilibriumSolver', defaultEquilibriumSolver, 'FLAG_RESULTS', false);
            defaultJumpConditionsSolver = combustiontoolbox.shockdetonation.JumpConditionsSolver('equilibriumSolver', defaultEquilibriumSolver, 'shockSolver', defaultShockSolver,'FLAG_RESULTS', false);
            defaultShockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical();
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultPlotConfig.plotProperties = {'K', 'R11', 'RTT', 'Ka', 'Kr'};

            % Parse inputs
            p = inputParser;
            addOptional(p, 'problemType', defaultProblemType, @(x) ischar(x) && any(strcmpi(x, {'compressible', 'vortical', 'acoustic', 'vortical_entropic'})));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver);
            addParameter(p, 'shockSolver', defaultShockSolver, @(x) isa(x, 'combustiontoolbox.shockdetonation.ShockSolver'));
            addParameter(p, 'jumpConditionsSolver', defaultJumpConditionsSolver, @(x) isa(x, 'combustiontoolbox.shockdetonation.JumpConditionsSolver'));
            addParameter(p, 'shockTurbulenceModel', defaultShockTurbulenceModel, @(x) isa(x, 'combustiontoolbox.shockturbulence.ShockTurbulenceModel'));
            addParameter(p, 'FLAG_TCHEM_FROZEN', obj.FLAG_TCHEM_FROZEN, @islogical);
            addParameter(p, 'FLAG_FROZEN', obj.FLAG_FROZEN, @islogical);
            addParameter(p, 'FLAG_INTERPOLATE', obj.FLAG_INTERPOLATE, @islogical);
            addParameter(p, 'FLAG_PAPER', obj.FLAG_PAPER, @islogical);
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set properties
            obj.problemType = p.Results.problemType;
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.shockSolver = p.Results.shockSolver;
            obj.jumpConditionsSolver = p.Results.jumpConditionsSolver;
            obj.shockTurbulenceModel = p.Results.shockTurbulenceModel;
            obj.FLAG_TCHEM_FROZEN = p.Results.FLAG_TCHEM_FROZEN;
            obj.FLAG_FROZEN = p.Results.FLAG_FROZEN;
            obj.FLAG_INTERPOLATE = p.Results.FLAG_INTERPOLATE;
            obj.FLAG_PAPER = p.Results.FLAG_PAPER;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;

            if ~sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.shockSolver.equilibriumSolver = obj.equilibriumSolver;
            else
                obj.equilibriumSolver.FLAG_TCHEM_FROZEN = obj.FLAG_TCHEM_FROZEN;
                obj.equilibriumSolver.FLAG_FROZEN = obj.FLAG_FROZEN;
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
            %     obj (EquilibriumSolver): Object of the class EquilibriumSolver
            
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
            %     obj (ShockSolver): ShockSolver object
            %     mixArray1 (Mixture): Array of Mixture objects (pre-shock state)
            %     mixArray2 (Mixture): Array of Mixture objects (post-shock state)
            %
            % Optional args:
            %     * mixArray_i (Mixture): Array of Mixture objects
            %
            % Example:
            %     * report(ShockTurbulenceSolver(), results);

            obj.plot(results);
        end

    end

end

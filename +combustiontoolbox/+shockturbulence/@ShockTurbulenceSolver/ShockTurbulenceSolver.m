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
    % These models are based on our previous theoretical works [1-3]
    % and have been extended to multi-component mixtures [4-7] using the
    % Combustion Toolbox [8-9].
    %
    % References:
    %     [1] Huete, C., Cuadra, A., Vera, M., Urzay, & J. (2021). Thermochemical
    %         effects on hypersonic shock waves interacting with weak turbulence.
    %         Physics of Fluids 33, 086111 (featured article). DOI: 10.1063/5.0059948.
    %
    %     [2] Huete, C., Velikovich, A. L., & Wouchuk, J. G. (2011). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional
    %         random isotropic density field. Physical Review E—Statistical, Nonlinear, and
    %         Soft Matter Physics, 83(5), 056320. DOI: 10.1103/PhysRevE.83.056320.
    %
    %     [3] Huete, C., Wouchuk, J. G., & Velikovich, A. L. (2012). Analytical linear theory
    %         for the interaction of a planar shock wave with a two-or three-dimensional random
    %         isotropic acoustic wave field. Physical Review E—Statistical, Nonlinear, and Soft
    %         Matter Physics, 85(2), 026312. DOI: 10.1063/5.0059948.
    %
    %     [4] Cuadra, A., Vera, M., Di Renzo, M., & Huete, C. (2023). Linear Theory
    %         of Hypersonic Shocks Interacting with Turbulence in Air. In 2023 AIAA
    %         SciTech Forum, National Harbor, USA. DOI: 10.2514/6.2023-0075.
    %
    %     [5] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %         and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %         Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %         Stanford University.
    %
    %     [6] Cuadra, A., Williams, C. T., Di Renzo, M., & Huete, C. The role of compressibility
    %         and vibrational-excitation in hypersonic shock–turbulence interactions.
    %         Journal of Fluid Mechanics (under review).
    %
    %     [7] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %         Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %         DOI: 10.1063/5.0255816.
    %
    %     [8] Cuadra, A., Huete, C., & Vera, M. (2026). Combustion Toolbox: An open-source
    %         thermochemical code for gas-and condensed-phase problems involving chemical equilibrium. 
    %         Computer Physics Communications 320, 110004. DOI:10.1016/j.cpc.2025.110004.
    %
    %     [9] Cuadra, A., Huete, C., Vera, M. (2022). Combustion Toolbox:
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
            defaultPlotConfig.plotProperties = {'K', 'R11', 'RTT', 'Ka', 'Kr', 'enstrophy', 'kolmogorovLengthRatio'};

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

            % Remove kolmogorovLengthRatio from plotConfig.plotProperties if problemType is 'acoustic'
            if strcmpi(obj.problemType, 'acoustic')
                obj.plotConfig.plotProperties = setdiff(obj.plotConfig.plotProperties, 'kolmogorovLengthRatio');
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

        function obj = set(obj, property, value, varargin)
            % Set properties of the ShockTurbulenceSolver object
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object with updated properties
            %
            % Example:
            %     set(ShockTurbulenceSolver(), 'problemType', 'compressible');
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function obj = setShockTurbulenceModel(obj, shockTurbulenceModel)
            % Set the ShockTurbulenceModel object or the name of the model
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     shockTurbulenceModel (ShockTurbulenceModel or char): ShockTurbulenceModel object or name of the model
            %
            % Returns:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object with updated ShockTurbulenceModel

            if isa(shockTurbulenceModel, 'combustiontoolbox.shockturbulence.ShockTurbulenceModel')
                obj.shockTurbulenceModel = shockTurbulenceModel;
                return
            end
                
            switch lower(shockTurbulenceModel)
                case 'vortical'
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVortical();
                case 'vortical_entropic'
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelVorticalEntropic();
                case 'acoustic'
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelAcoustic();
                case 'compressible'
                    obj.shockTurbulenceModel = combustiontoolbox.shockturbulence.ShockTurbulenceModelCompressible();
                otherwise
                    error('Unknown shock turbulence model: %s', shockTurbulenceModel);
            end
            
        end

        function results = solve(obj, mixArray1, varargin)
            % Solve array of shock waves problems
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     mixArray1 (Mixture): Initial Mixture objects
            %
            % Returns:
            %     results (struct): Struct with results from LIA
            %
            % Example:
            %     * results = solve(ShockTurbulenceSolver(), mixArray1);
            
            % Definitions
            eta = [mixArray1.eta];
            etaVorticity = [mixArray1.etaVorticity];
            chi = [mixArray1.chi];
            viscosityModel = obj.shockTurbulenceModel.viscosityModel;

            % Timer
 	        obj.time = tic;

            % Compute jump conditions
            [jumpConditions, mixArray1, mixArray2] = obj.getJumpConditions(mixArray1, varargin{:});
            
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
                    results = obj.shockTurbulenceModel.getAverages(Rratio, M2, Gammas2, Gammas1, Gammas3, beta, eta, chi, etaVorticity);
            end

            % Add jump conditions
            names = fieldnames(jumpConditions);
            for k = 1:length(names)
                results.(names{k}) = jumpConditions.(names{k});
            end

            % Compute Kolmogorov length scale ratio across the shock
            if ~strcmpi(obj.problemType, 'acoustic')
                results.kolmogorovLengthRatio = obj.getKolmogorovLength(results, 'viscosityModel', viscosityModel);
            end

            % Assign mixArray1 and mixArray2 to results
            results.mixArray1 = mixArray1;
            results.mixArray2 = mixArray2;

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
            ax1 = plotProperties(repmat({results.mixArray1(1).rangeName}, 1, numPlotProperties), results.mixArray1, obj.plotConfig.plotProperties, results, 'config', obj.plotConfig);
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

        function kolmogorovLengthRatio = getKolmogorovLength(obj, results, varargin)
            % Estimate Kolmogorov length scale ratio across the shock
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     results (struct): Struct with results from LIA
            %
            % Optional name-value pairs:
            %     * viscosityModel (char): Viscosity model to compute dynamic viscosity ratio across the shock ('powerlaw' or 'sutherland')
            %
            % Returns:
            %     kolmogorovLengthRatio (float): Kolmogorov length scale ratio
            %
            % Example:
            %     kolmogorovLengthRatio = getKolmogorovLength(ShockTurbulenceSolver(), results);
            %
            % Note: The calculation of the dynamic viscosity ratio is based on temporal functions and will be overriden with a specific TransportProperties class in future releases.

            % Check STI model
            if strcmpi(obj.problemType, 'acoustic')
                error('Kolmogorov length scale ratio can only be computed for vortical, vortical-entropic or compressible disturbances');
            end

            % Definitions
            T1 = results.T1;
            Rratio = results.Rratio;
            TRatio = results.Tratio;
            enstrophyRatio = results.enstrophy;
            defaultViscosityModel = 'powerlaw';

            % Select viscosity model
            p = inputParser;
            addParameter(p, 'viscosityModel', defaultViscosityModel, @(x) ischar(x) && any(strcmpi(x, {'powerlaw', 'sutherland'})));
            parse(p, varargin{:});
            viscosityModel = p.Results.viscosityModel;
            
            % Compute dynamic viscosity ratio across the shock
            switch lower(viscosityModel)
                case 'powerlaw'
                    muRatio = obj.getDynamicViscosityPowerLawRatio(T1, TRatio .* T1);
                case 'sutherland'
                    muRatio = obj.getDynamicViscositySutherlandRatio(T1, TRatio .* T1);
                otherwise
                    error('Unknown viscosity model: %s', viscosityModel);
            end

            % Kolmogorov length scale ratio 
            kolmogorovLengthRatio = Rratio.^(-1/2) .* muRatio.^(1/2) .* enstrophyRatio.^(-1/4);
        end

    end

    methods (Access = private)

        function [jumpConditions, mixArray1, mixArray2] = getJumpConditions(obj, mixArray1, varargin)
            % Get jump conditions across the shock for an array of Mixture objects
            %
            % Args:
            %     obj (ShockTurbulenceSolver): ShockTurbulenceSolver object
            %     mixArray1 (Mixture): Initial Mixture objects
            %
            % Returns:
            %     jumpConditions (struct): Struct with jump conditions across the shock
            %     mixArray1 (Mixture): Pre-shock Mixture objects
            %     mixArray2 (Mixture): Post-shock Mixture objects

            % Check that mixArray is not evaluated at the same upstream conditions (T, p, and mach)
            if length(unique([mixArray1.T])) ~= 1 || length(unique([mixArray1.p])) ~= 1 || length(unique([mixArray1.mach])) ~= 1
                
                [jumpConditions, mixArray1, mixArray2] = obj.jumpConditionsSolver.solve(mixArray1, varargin{:});
                return
            end

            [jumpConditions, ~, mixArray2] = obj.jumpConditionsSolver.solve(mixArray1(1), varargin{:});

            % Copy all fields inside jumpConditions to match size of mixArray
            jumpConditionsFields = fieldnames(jumpConditions);
            for i = length(jumpConditionsFields):-1:1
                name = jumpConditionsFields{i};
                jumpConditions.(name) = repmat(jumpConditions.(name), size(mixArray1));
            end
            
            % Copy all entries to match size of mixArray
            mixArray2 = repmat(mixArray2, size(mixArray1));
        end

    end

    methods (Access = private, Static)
    
        function muRatio = getDynamicViscosityPowerLawRatio(T1, T2)
            % Get dynamic viscosity ratio using a power-law model
            %
            % Args:
            %     T1 (float): Pre-shock temperature [K]
            %     T2 (float): Post-shock temperature [K]
            %
            % Returns:
            %     muRatio (float): Dynamic viscosity ratio mu2/mu1
            %
            % Note: This is a temporal function and will be overriden with a specific TransportProperties class in future releases.

            % Compute dynamic viscosity ratio using power-law
            muRatio = (T2 ./ T1).^(3/4);
        end

        function muRatio = getDynamicViscositySutherlandRatio(T1, T2)
            % Get dynamic viscosity ratio using Sutherland's law model
            %
            % Args:
            %     T1 (float): Pre-shock temperature [K]
            %     T2 (float): Post-shock temperature [K]
            %
            % Returns:
            %     muRatio (float): Dynamic viscosity ratio mu2/mu1
            %
            % Note: This is a temporal function and will be overriden with a specific TransportProperties class in future releases.

            % Definitions
            S = 110.4; % Sutherland's constant [K]

            % Compute dynamic viscosity ratio using Sutherland's law
            muRatio = (T2 ./ T1).^(3/2) .* (T1 + S) ./ (T2 + S);
        end

    end

end

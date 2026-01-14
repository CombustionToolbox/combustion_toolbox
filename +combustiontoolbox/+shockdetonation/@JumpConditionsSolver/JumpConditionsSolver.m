classdef JumpConditionsSolver < handle
    % The :mat:func:`JumpConditionsSolver` class computes the Rankine–Hugoniot 
    % jump conditions and associated amplification parameters for gas mixtures 
    % across a shock wave. It supports calorically perfect, thermally perfect, and 
    % calorically imperfect models. In addition to post-shock thermodynamic properties, 
    % the class evaluates the dimensionless slopes of the Hugoniot curve—``Gammas2`` 
    % (:math:`\Gamma`), ``Gammas1`` (:math:`\Gamma_\rho`), and ``Gammas3`` (:math:`\Gamma_p`)—which 
    % quantify the local sensitivity of the post-shock density to changes in downstream 
    % pressure (:math:`\Gamma`), and the influence of upstream density and pressure 
    % perturbations on the post-shock state (:math:`\Gamma_\rho`, :math:`\Gamma_p`), 
    % as formalized in Ref. [1, 2].
    %
    % The :mat:func:`JumpConditionsSolver` object can be initialized as follows: ::
    %
    %       solver = JumpConditionsSolver('equilibriumSolver', eqSolver, ...)
    %
    % Optional name-value arguments can be specified to customize tolerances, 
    % interpolation settings, plotting options, and solver behavior. Default 
    % objects for :mat:func:`EquilibriumSolver` and :mat:func:`ShockSolver` 
    % are created automatically if not provided.
    %
    % Flags:
    %     * ``FLAG_FAST``: Use previous state as initial guess for equilibrium solver
    %     * ``FLAG_PLOT``: Plot results
    %     * ``FLAG_PAPER``: Apply the normalization and notation conventions described in Cuadra2025a
    %     * ``FLAG_RESULTS``: Print computed values to console
    %     * ``FLAG_TIME``: Print elapsed computation time
    %     * ``FLAG_REPORT``: Generate predefined plots for diagnostic analysis
    %
    % See also: :mat:func:`ShockSolver`, :mat:func:`EquilibriumSolver`, 
    % :mat:func:`solve`, :mat:func:`report`, :mat:func:`plot`, :mat:func:`plotGammas`
    %
    % References:
    %
    %      [1] Cuadra, A., Di Renzo, M., Hoste, J. J. O., Williams, C. T., Vera, M., & Huete, C. (2025).
    %          Review of shock-turbulence interaction with a focus on hypersonic flow. Physics of Fluids, 37(4).
    %          DOI: 10.1063/5.0255816.
    %
    %      [2] Cuadra, A., Williams, C. T., Di Renzo, M. & Huete, C. (2024). Compressibility
    %          and vibrational-excitation effects in hypersonic shock-turbulence interaction.
    %          Tech. Rep. Summer Program Proceedings, Center for Turbulence Research,
    %          Stanford University.

    properties
        equilibriumSolver              % EquilibriumSolver object
        shockSolver                    % ShockSolver object
        tolGammas1 = 1e-4;             % Tolerance for the calculation of the dimensionless RH slope Gammas1
        tolGammas2 = 1e-4;             % Tolerance for the calculation of the dimensionless RH slope Gammas2
        tolGammas3 = 1e-4;             % Tolerance for the calculation of the dimensionless RH slope Gammas3
        FLAG_PLOT = false              % Flag to plot results
        FLAG_PAPER = false             % Flag to apply the normalization and notation conventions described in Cuadra2025a
        FLAG_RESULTS = true            % Flag to print results
        FLAG_TIME = true               % Flag to print elapsed time
        FLAG_REPORT = false            % Flag to print predefined plots
        time                           % Elapsed time
        plotConfig
    end

    properties (Access = private)
        mix                            % Mixture object for initialization
        temperature                    % Temperature array of the array of Mixture objects
        pressure                       % Pressure array of the array of Mixture objects
        mach                           % Mach number array of the array of Mixture objects
        dp2dM1                         % Derivative of p2 with respect to M1
        drho2dM1                       % Derivative of rho2 with respect to M1
    end

    properties (Dependent, Access = private)
        numCases
        FLAG_SCALAR
    end

    methods 

        function obj = JumpConditionsSolver(varargin)
            % Constructor
            % 
            % Optional Args:
            %     * param1 (type): Description of param1
            %     * param2 (type): Description of param2
            %
            % Returns:
            %     obj (JumpConditionsSolver): Instance of the solver
            %
            % Example:
            %     obj = JumpConditionsSolver('param1', value1, 'param2', value2, ...);

            % Default values
            defaultCaloricGasModel = combustiontoolbox.core.CaloricGasModel.imperfect;
            defaultEquilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver();
            defaultShockSolver = combustiontoolbox.shockdetonation.ShockSolver('equilibriumSolver', defaultEquilibriumSolver);
            defaultFLAG_TCHEM_FROZEN = false;
            defaultFLAG_FROZEN = false;
            defaultFLAG_FAST = true;
            defaultPlotConfig = combustiontoolbox.utils.display.PlotConfig();
            defaultPlotConfig.innerposition = [0.15 0.15 0.35 0.45];
            defaultPlotConfig.outerposition = [0.15 0.15 0.35 0.45];
            defaultPlotConfig.plotProperties = {'Rratio', 'Pratio', 'Tratio', 'M2'};

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'caloricGasModel', defaultCaloricGasModel, @(x) isa(x, 'combustiontoolbox.core.CaloricGasModel'));
            addParameter(p, 'equilibriumSolver', defaultEquilibriumSolver, @(x) isa(x, 'combustiontoolbox.equilibrium.EquilibriumSolver'));
            addParameter(p, 'shockSolver', defaultShockSolver, @(x) isa(x, 'combustiontoolbox.shockdetonation.ShockSolver'));
            addParameter(p, 'tolGammas1', obj.tolGammas1, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'tolGammas2', obj.tolGammas2, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'tolGammas3', obj.tolGammas3, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'FLAG_FAST', defaultFLAG_FAST, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_PLOT', obj.FLAG_PLOT, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_PAPER', obj.FLAG_PAPER, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_RESULTS', obj.FLAG_RESULTS, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_TIME', obj.FLAG_TIME, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_REPORT', obj.FLAG_REPORT, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_TCHEM_FROZEN', defaultFLAG_TCHEM_FROZEN, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'FLAG_FROZEN', defaultFLAG_FROZEN, @(x) islogical(x) && isscalar(x));
            addParameter(p, 'plotConfig', defaultPlotConfig, @(x) isa(x, 'combustiontoolbox.utils.display.PlotConfig'));
            parse(p, varargin{:});

            % Set parameters
            obj.equilibriumSolver = p.Results.equilibriumSolver;
            obj.shockSolver = p.Results.shockSolver;
            obj.tolGammas1 = p.Results.tolGammas1;
            obj.tolGammas2 = p.Results.tolGammas2;
            obj.tolGammas3 = p.Results.tolGammas3;
            obj.FLAG_PLOT = p.Results.FLAG_PLOT;
            obj.FLAG_PAPER = p.Results.FLAG_PAPER;
            obj.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            obj.FLAG_TIME = p.Results.FLAG_TIME;
            obj.FLAG_REPORT = p.Results.FLAG_REPORT;
            obj.plotConfig = p.Results.plotConfig;

            if sum(contains(p.UsingDefaults, 'equilibriumSolver'))
                obj.equilibriumSolver.FLAG_FAST = p.Results.FLAG_FAST;
                obj.equilibriumSolver.caloricGasModel = p.Results.caloricGasModel;
            end

            if sum(contains(p.UsingDefaults, 'shockSolver'))
                obj.shockSolver.FLAG_RESULTS = p.Results.FLAG_RESULTS;
            end

            % Miscellaneous
            obj.equilibriumSolver.FLAG_RESULTS = false;
            obj.equilibriumSolver.FLAG_TIME = false;
            obj.equilibriumSolver.FLAG_CACHE = false;
            obj.shockSolver.FLAG_TIME = false;

            % Display warning if deprecated flags are used
            if ~ismember('FLAG_TCHEM_FROZEN', p.UsingDefaults) || ~ismember('FLAG_FROZEN', p.UsingDefaults)
                warning(['The flags ''FLAG_TCHEM_FROZEN'' and ''FLAG_FROZEN'' are deprecated. ', ...
                         'Please use the ''caloricGasModel'' parameter with values from the CaloricGasModel enumeration instead.']);
            
                obj.equilibriumSolver.caloricGasModel = obj.equilibriumSolver.caloricGasModel.fromFlag(p.Results.FLAG_TCHEM_FROZEN, p.Results.FLAG_FROZEN);
            end

            % Assign equilibriumSolver to shockSolver
            obj.shockSolver.equilibriumSolver = obj.equilibriumSolver;
        end

        function value = get.numCases(obj)
            % Get number of cases
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %
            % Returns:
            %     value (int): Number of cases
            
            value = length(obj.mach);
        end

        function value = get.FLAG_SCALAR(obj)
            % Get flag indicating if the problem is scalar
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %
            % Returns:
            %     value (bool): Flag indicating if the problem is scalar
            
            value = isscalar(obj.mach);
        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the JumpConditionsSolver object
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object with updated properties
            %
            % Examples:
            %     * set(JumpConditionsSolver(), 'tolGammas1', 1e-5);
            %     * set(JumpConditionsSolver(), 'FLAG_PAPER', true);
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

        function [jumpConditions, mixArray1, mixArray2] = solve(obj, mixArray1, varargin)
            % Solve the jump conditions for the given mixture
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     mixArray (Mixture): Array of Mixture objects
            %
            % Returns:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object with updated properties
            %
            % Example:
            %     jumpConditions = solve(JumpConditionsSolver(), mixArray);
            
            % Definitions
            obj.temperature = [mixArray1.T];
            obj.pressure = [mixArray1.p];
            obj.mach = [mixArray1.mach];
            obj.mix = copy(mixArray1(1)); obj.mix.u = [];
            
            % Timer
            obj.time = tic;
            
            % Compute jump conditions
            [jumpConditions, mixArray1, mixArray2] = getJumpData(obj, mixArray1);

            % Solve Gammas (calorically perfect gas)
            if obj.equilibriumSolver.caloricGasModel.isPerfect()
                [jumpConditions.Gammas1, jumpConditions.Gammas2, jumpConditions.Gammas3] = obj.getGammasPerfect(jumpConditions.gamma1, jumpConditions.M1);

                % Timer
                obj.time = toc(obj.time);

                % Print elapsed time
                printTime(obj);
                return
            end

            % Calculate Gammas1, Gammas2, and Gammas3
            jumpConditions.Gammas2 = getGammas2(obj, jumpConditions);
            jumpConditions.Gammas1 = getGammas1(obj, jumpConditions);
            jumpConditions.Gammas3 = getGammas3(obj, jumpConditions);

            % Timer
            obj.time = toc(obj.time);

            % Print elapsed time
            printTime(obj);

            % Postprocess all the results with predefined plots
            if obj.FLAG_REPORT
                report(obj, jumpConditions);
            end

        end

        function report(obj, jumpConditions)
            % Generate report for the jump conditions
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Example:
            %     report(JumpConditionsSolver(), jumpConditions);

            % Plot results
            plot(obj, jumpConditions);
            plotGammas(obj, jumpConditions);
        end

        function ax1 = plot(obj, jumpConditions)
            % Plot results
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Example:
            %     * plot(JumpConditionsSolver(), results);
            
            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            numPlotProperties = obj.plotConfig.numPlotProperties;

            % Check if is a scalar value
            if isscalar(jumpConditions.M1) || numPlotProperties < 1
                ax1 = [];
                return
            end
            
            % Plot properties
            ax1 = plotProperties(repmat({'M1'}, 1, numPlotProperties), jumpConditions, obj.plotConfig.plotProperties, jumpConditions, 'config', obj.plotConfig);
        end

        function ax = plotGammas(obj, jumpConditions)
            % Plot dimensionless slopes of the Hugoniot curve
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Returns:
            %     ax (axes): Axes handle for the plot
            %
            % Example:
            %     plotGammas(JumpConditionsSolver(), jumpConditions);
        
            % Import packages
            import combustiontoolbox.utils.display.*
            
            % Definitions
            M1 = jumpConditions.M1;
            Gammas1 = jumpConditions.Gammas1;
            Gammas2 = jumpConditions.Gammas2;
            Gammas3 = jumpConditions.Gammas3;

            if obj.FLAG_SCALAR
                return
            end

            % Get calorically perfect gas values
            [Gammas1_perfect, Gammas2_perfect, Gammas3_perfect] = obj.getGammasPerfect(jumpConditions.gamma1, M1);

            % Set figure
            ax = setFigure(obj.plotConfig);
            
            % Plot figure
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', -Gammas1, 'ax', ax, 'config', obj.plotConfig, 'color', obj.plotConfig.gray, 'linestyle', '-');
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', Gammas2, 'ax', ax, 'config', obj.plotConfig, 'color', obj.plotConfig.blue, 'linestyle', '-');
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', Gammas3, 'ax', ax, 'config', obj.plotConfig, 'color', [0 0 0], 'linestyle', '-');
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', -Gammas1_perfect, 'ax', ax, 'config', obj.plotConfig, 'color', obj.plotConfig.gray, 'linestyle', '--');
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', Gammas2_perfect, 'ax', ax, 'config', obj.plotConfig, 'color', obj.plotConfig.blue, 'linestyle', '--');
            plotFigure(' $ Pre-shock Mach number, $\mathcal{M}_1', M1, ' $ Dimensionless Hugoniot-slope, $\Gamma_i', Gammas3_perfect, 'ax', ax, 'config', obj.plotConfig, 'color', [0 0 0], 'linestyle', '--');
            
            if obj.FLAG_PAPER
                legend(ax, {'$-\Gamma_\rho$', '$\Gamma$', '$\Gamma_p$', '$-\Gamma_\rho^{\rm perfect}$', '$\Gamma^{\rm perfect}$', '$\Gamma_p^{\rm perfect}$'}, 'Interpreter', 'latex', 'Location', 'best');
                return
            end

            legend(ax, {'$-\Gamma_1$', '$\Gamma_2$', '$\Gamma_3$', '$-\Gamma_1^{\rm perfect}$', '$\Gamma_2^{\rm perfect}$', '$\Gamma_3^{\rm perfect}$'}, 'Interpreter', 'latex', 'Location', 'best');
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

    methods (Access = private)

        function [jumpConditions, mixArray1, mixArray2] = getJumpData(obj, mixArray1)
            % Get jump conditions
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     mixArray1 (Mixture): Mixture object for the pre-shock state
            %
            % Returns:
            %     Tuple containing:
            %
            %     * jumpConditions (struct): Structure containing the jump conditions
            %     * mixArray1 (Mixture): Mixture object for the pre-shock state
            %     * mixArray2 (Mixture): Mixture object for the post-shock state

            % Print status
            fprintf('Solving jump conditions... ');

            % Initialize Gammas
            Gammas1 = [];
            Gammas2 = [];
            Gammas3 = [];

            % Solve jump conditions
            [mixArray1, mixArray2] = obj.shockSolver.solveArray(mixArray1);

            % Get properties
            rho1 = [mixArray1.rho]; rho2 = [mixArray2.rho];
            p1 = [mixArray1.p]; p2 = [mixArray2.p];
            T1 = [mixArray1.T]; T2 = [mixArray2.T];
            M1 = [mixArray1.mach]; M2 = [mixArray2.mach];
            a1 = [mixArray1.sound]; a2 = [mixArray2.sound];
            u1 = [mixArray1.u]; u2 = [mixArray2.uShock];
            m1 = [mixArray1.mi]; m2 = [mixArray2.mi];
            gamma1 = [mixArray1.gamma_s]; gamma2 = [mixArray2.gamma_s];
            Rg1 = ([mixArray1.cp] - [mixArray1.cv]) ./ m1; 
            Rg2 = ([mixArray2.cp] - [mixArray2.cv]) ./ m2;

            % Compute jump conditions
            Rratio = rho2 ./ rho1;
            Pratio = p2 ./ p1;
            Tratio = T2 ./ T1;
            beta = a2 ./ a1;

            % Set jumpConditions
            jumpConditions.Rratio = Rratio;
            jumpConditions.Pratio = Pratio;
            jumpConditions.Tratio = Tratio;
            jumpConditions.M1 = M1;
            jumpConditions.M2 = M2;
            jumpConditions.Gammas1 = Gammas1;
            jumpConditions.Gammas2 = Gammas2;
            jumpConditions.Gammas3 = Gammas3;
            jumpConditions.beta = beta;
            jumpConditions.rho1 = rho1; jumpConditions.rho2 = rho2;
            jumpConditions.p1 = p1; jumpConditions.p2 = p2;
            jumpConditions.T1 = T1; jumpConditions.T2 = T2;
            jumpConditions.a1 = a1; jumpConditions.a2 = a2;
            jumpConditions.u1 = u1; jumpConditions.u2 = u2;
            jumpConditions.m1 = m1; jumpConditions.m2 = m2;
            jumpConditions.gamma1 = gamma1; jumpConditions.gamma2 = gamma2;
            jumpConditions.Rg1 = Rg1; jumpConditions.Rg2 = Rg2;

            % Print status
            fprintf('OK!\n');
        end

        function [rho1, rho2, p2] = solveShock(obj, T1, p1)
            % Solve shock problem for the given temperature and pressure
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     T1 (float): Temperature of the mixture
            %     p1 (float): Pressure of the mixture
            %
            % Returns:
            %     Tuple containing:
            %
            %     * rho1 (float): Density of the mixture in the pre-shock state
            %     * rho2 (float): Density of the mixture in the post-shock state
            %     * p2 (float): Pressure of the mixture in the post-shock state

            % Definitions
            M1 = obj.mach;

            % Solve base solution
            mixArray1 = setProperties(obj.mix, 'temperature', T1, 'pressure', p1, 'M1', M1);
        
            % Solve problem
            [mixArray1, mixArray2] = obj.shockSolver.solveArray(mixArray1);

            % Get properties
            rho1 = [mixArray1.rho];
            rho2 = [mixArray2.rho];
            p2 = [mixArray2.p];
        
        end

        function [Gammas1, Gammas2, Gammas3] = getGammasPerfect(obj, gamma, mach)
            % Get the dimensionless slope of the Hugoniot curve for a calorically perfect gas
            %
            % Args:
            %     gamma (float): Adiabatic index
            %     mach (float): mach number
            %
            % Returns:
            %     Tuple containing:
            %
            %     * Gammas1 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho2, p1)
            %     * Gammas2 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, p1)
            %     * Gammas3 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, rho2)
            
            % Definitions
            R = obj.getDensityRatioPerfect(gamma, mach);
            P = obj.getPressureRatioPerfect(gamma, mach);
            
            % Compute Gammas
            Gammas2 = 1 ./ mach.^2;
            Gammas1 = - R .* Gammas2;
            Gammas3 =  Gammas2;

            if ~obj.FLAG_PAPER
                return
            end

            % Compute Gammas using Cuadra2024b and Cuadra2025a formulae
            Gammas1 = Gammas2 ./ Gammas1;
            Gammas3 = Gammas2.^2 ./ Gammas3 .* P;
        end
        
        function Gammas1 = getGammas1(obj, jumpConditions)
            % Get the dimensionless slope of the Hugoniot curve Gammas1
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Returns:
            %     Gammas1 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho2, p1)
            
            % Import packages
            import combustiontoolbox.common.Units

            % Print status
            fprintf('Solving Gammas1... ');

            % Definitions
            u1 = jumpConditions.u1;
            rho1 = jumpConditions.rho1;
            p2 = jumpConditions.p2;
            rho2 = jumpConditions.rho2;
            Gammas2 = jumpConditions.Gammas2;

            % Perturb temperature
            T1_Gammas1 = obj.temperature .* (1 + obj.tolGammas1);

            % Calculate jump conditions for the perturbed temperature
            [rho1_Gammas1, rho2_Gammas1, p2_Gammas1] = obj.solveShock(T1_Gammas1, obj.pressure);

            % Get derivates of p2 and rho2 with respect to M1
            dp2dM1 = obj.dp2dM1;
            drho2dM1 = obj.drho2dM1;

            % Compute derivatives with respect to T1
            dT1 = (T1_Gammas1 - obj.temperature);
            dp2dT1 = Units.bar2Pa * (p2_Gammas1 - p2) ./ dT1;
            drho2dT1 = (rho2_Gammas1 - rho2) ./ dT1;

            % Compute dp2/dT1 at constant rho2 using the chain rule
            dp2dT1_rho2 = dp2dT1 - dp2dM1 .* (drho2dT1 ./ drho2dM1);

            % Compute Gammas1
            Gammas1 = u1.^2 .* (rho1_Gammas1 - rho1) ./ (dp2dT1_rho2 .* dT1);

            if obj.FLAG_PAPER
                Gammas1 = Gammas2 ./ Gammas1;
            end

            % Print status
            fprintf('OK!\n');
        end

        function Gammas2 = getGammas2(obj, jumpConditions)
            % Get the dimensionless slope of the Hugoniot curve Gammas2
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Returns:
            %     Gammas2 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, p1)
            
            % Import packages
            import combustiontoolbox.common.Units

            % Print status
            fprintf('Solving Gammas2... ');

            % Definitions
            T1 = obj.temperature;
            p1 = obj.pressure;
            M1 = obj.mach;

            % Perturb Mach number
            M1_Gammas2 = M1 .* (1 + obj.tolGammas2);

            % Define new mixture
            mixArray1 = setProperties(obj.mix, 'temperature', T1, 'pressure', p1, 'M1', M1_Gammas2);

            % Calculate jump conditions for the perturbed temperature
            [~, mixArray2] = obj.shockSolver.solveArray(mixArray1);

            % Get properties
            rho2 = obj.reshapeVector(jumpConditions.rho2, [mixArray2.rho]);
            p2 = obj.reshapeVector(jumpConditions.p2, [mixArray2.p]);
            u2 = obj.reshapeVector(jumpConditions.u2, [mixArray2.uShock]);
            
            % Compute Gammas2
            p2_Pa = p2 * Units.bar2Pa;
            indexBase = 1:2:numel(p2_Pa)-1;
            indexPert = indexBase + 1;
            Gammas2 = u2(indexBase).^2 .* (rho2(indexPert) - rho2(indexBase)) ./ (p2_Pa(indexPert) - p2_Pa(indexBase));

            % Compute derivatives of p2 and rho2 with respect to M1
            obj.dp2dM1 = (p2_Pa(indexPert) - p2_Pa(indexBase)) ./ (M1_Gammas2 - M1);
            obj.drho2dM1 = (rho2(indexPert) - rho2(indexBase)) ./ (M1_Gammas2 - M1);

            % Print status
            fprintf('OK!\n');
        end

        function Gammas3 = getGammas3(obj, jumpConditions)
            % Get the dimensionless slope of the Hugoniot curve Gammas3
            %
            % Args:
            %     obj (JumpConditionsSolver): JumpConditionsSolver object
            %     jumpConditions (struct): Structure containing the jump conditions
            %
            % Returns:
            %     Gammas3 (float): Dimensionless slope of the Hugoniot curve (partial derivative at constant rho1, rho2)
            
            % Import packages
            import combustiontoolbox.common.Units

            % Print status
            fprintf('Solving Gammas3... ');

            % Definitions
            Rg1 = jumpConditions.Rg1;
            rho1 = jumpConditions.rho1;
            rho2 = jumpConditions.rho2;
            p1 = jumpConditions.p1;
            p2 = jumpConditions.p2;
            Rratio = jumpConditions.Rratio;
            beta = jumpConditions.beta;
            Gammas2 = jumpConditions.Gammas2;

            % Perturb pressure
            p1_Gammas3 = obj.pressure .* (1 + obj.tolGammas3);

            % Correct temperature to remain the density constant
            T1_Gammas3 = p1_Gammas3 .* Units.bar2Pa ./ (rho1 .* Rg1);

            % Calculate jump conditions for the perturbed temperature
            [~, rho2_Gammas3, p2_Gammas3] = obj.solveShock(T1_Gammas3, p1_Gammas3);

            % Get derivates of p2 and rho2 with respect to M1
            dp2dM1 = obj.dp2dM1;
            drho2dM1 = obj.drho2dM1;

            % Compute derivatives with respect to p1
            dp2dp1 = (p2_Gammas3 - p2) ./ (p1_Gammas3 - p1);
            drho2dp1 = (rho2_Gammas3 - rho2) ./ (Units.bar2Pa * (p1_Gammas3 - p1));

            % Compute dp2/dp1 at constant rho2 using the chain rule
            dp2dp1_rho2 = dp2dp1 - dp2dM1 .* (drho2dp1 ./ drho2dM1);

            % Compute Gammas3
            Gammas3 = Gammas2 .* Rratio .* beta.^2 .* dp2dp1_rho2.^(-1);

            if obj.FLAG_PAPER
                Gammas3 = Gammas2 ./ Gammas3 .* (Gammas2 .* Rratio .* beta.^2);
            end

            % Print status
            fprintf('OK!\n');
        end


    end

    methods (Static, Access = public)

        function Rratio = getDensityRatioPerfect(gamma, mach)
            % Get the density ratio for a calorically perfect gas
            %
            % Args:
            %     gamma (float): Adiabatic index
            %     mach (float): mach number
            %
            % Returns:
            %     Rratio (float): Density ratio
            
            Rratio = ( (gamma + 1) .* mach.^2 ) ./ ( gamma .* mach.^2 + 1 - (mach.^2 - 1) );
        end

        function Pratio = getPressureRatioPerfect(gamma, mach)
            % Get the pressure ratio for a calorically perfect gas
            %
            % Args:
            %     gamma (float): Adiabatic index
            %     mach (float): mach number
            %
            % Returns:
            %     Pratio (float): Pressure ratio
            
            % Definitions
            Rratio = combustiontoolbox.shockdetonation.JumpConditionsSolver.getDensityRatioPerfect(gamma, mach);

            % Compute pressure ratio
            Pratio = 1 + gamma .* mach.^2 .* (1 - 1 ./ Rratio);
        end

    end

    methods (Static, Access = private)

        function vector = reshapeVector(x, y, varargin)
            % Reshape vectors to be column vectors of the form
            %
            %   [x(1), y(1), x(2), y(2), ...]
            %
            % Args:
            %     x (float): First vector
            %     y (float): Second vector
            %
            % Optional Args:
            %     * varargin (float): Additional vectors

            % Definitions
            vector = zeros(length(x), nargin);
            vector(:, 1) = x;
            vector(:, 2) = y;

            % Reshape vectors
            for i = 3:nargin
                vector(:, i) = varargin{i};
            end

            vector = reshape(vector.', 1, []);
        end

    end

end

classdef AppSolver < handle
    % Dispatches GUI study requests to CT solvers and packages results
    %
    % Attributes:
    %     session (AppSession): Long-lived GUI services
    %     problemBuilder (AppProblemBuilder): Builder for solver problem setup
    %     problemCatalog (AppProblemCatalog): Problem metadata and output layouts
    %
    % Example:
    %     appSolver = combustiontoolbox.gui.AppSolver(session);

    properties (Access = private)
        session
        problemBuilder
        problemCatalog
    end

    methods
        function obj = AppSolver(session, varargin)
            % AppSolver constructor
            %
            % Args:
            %     session (AppSession): Long-lived GUI services
            %
            % Optional Args:
            %     problemCatalog (AppProblemCatalog): Problem metadata catalog
            %
            % Returns:
            %     obj (AppSolver): Initialized GUI solver orchestrator
            if nargin == 0
                return
            end

            obj.session = session;
            obj.problemBuilder = combustiontoolbox.gui.AppProblemBuilder(session);
            obj.problemCatalog = combustiontoolbox.gui.AppProblemCatalog();

            if nargin > 1 && isa(varargin{1}, 'combustiontoolbox.gui.AppProblemCatalog')
                obj.problemCatalog = varargin{1};
            end
        end

        function problem = buildProblem(obj, input, varargin)
            % Build the solver-ready problem setup
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Optional Args:
            %     mixture (Mixture): Current mixture state from the GUI
            %
            % Returns:
            %     problem (struct): Solver-ready problem setup
            problem = obj.problemBuilder.build(input, varargin{:});
        end

        function script = exportScript(obj, setup)
            % Export canonical setup data to a MATLAB script
            %
            % Args:
            %     setup (struct): Canonical problem setup
            %
            % Returns:
            %     script (char): Runnable MATLAB script text
            script = obj.problemBuilder.exportScript(setup);
        end

        function solution = solve(obj, input, varargin)
            % Solve a GUI study request
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Optional Args:
            %     mixture (Mixture): Current mixture state from the GUI
            %
            % Returns:
            %     solution (struct): Problem setup, solver outputs, and GUI results
            problem = obj.buildProblem(input, varargin{:});
            setup = problem.setup;
            definition = obj.problemCatalog.get(setup.requestedProblemType);
            route = definition.solverRoute;
            layout = definition.resultLayout;
            solver = obj.solverFromRoute(route);
            solverProblemType = obj.solverProblemType(problem.problemType, route);

            obj.configureSolver(solver, solverProblemType, setup.options.printResults, definition.plotConfig, route);
            solverOutputs = obj.solveWithProblemOptions(solver, route, problem, layout);
            reportOutputs = obj.selectReportOutputs(route, solverOutputs);
            obj.configureReport(solver, setup.options);

            resultStore = combustiontoolbox.gui.AppResults();
            results = resultStore.buildFromSolverOutputs(problem.problemType, solverOutputs, layout, problem.problemType);
            results = obj.attachSetupMetadata(results, setup);

            solution = struct();
            solution.problem = problem;
            solution.solver = solver;
            solution.solverProblemType = solverProblemType;
            solution.outputVariant = problem.problemType;
            solution.resultLayout = layout;
            solution.solverOutputs = solverOutputs;
            solution.reportOutputs = reportOutputs;
            solution.results = results;
        end
    end

    methods (Access = private)
        function solver = solverFromRoute(obj, route)
            % Return the solver object associated with a catalog route
            %
            % Args:
            %     route (struct): Solver route metadata from AppProblemCatalog
            %
            % Returns:
            %     solver (handle): Solver instance stored in the GUI session
            if isempty(obj.session) || ~isprop(obj.session, route.solver) || isempty(obj.session.(route.solver))
                error('AppSolver:MissingSolver', ...
                    'The session does not provide the solver "%s".', route.solver);
            end

            solver = obj.session.(route.solver);
        end

        function configureSolver(obj, solver, solverProblemType, printResults, plotConfig, route)
            % Configure a solver before executing the current problem
            %
            % Args:
            %     solver (handle): Solver object used for the current problem
            %     solverProblemType (char): Solver-level problem type identifier
            %     printResults (logical): Flag to print solver reports
            %     plotConfig (struct): Plot configuration metadata from the catalog
            %     route (struct): Solver route metadata from AppProblemCatalog
            solver.problemType = solverProblemType;

            if isprop(solver, 'FLAG_RESULTS')
                solver.FLAG_RESULTS = printResults;
            end

            if strcmp(route.solver, 'shockTurbulenceSolver')
                obj.configureShockTurbulenceSolver(solver, solverProblemType);
            end

            if isprop(solver, 'plotConfig') && ~isempty(solver.plotConfig)
                solver.plotConfig.plotProperties = plotConfig.properties;
                solver.plotConfig.plotPropertiesBasis = plotConfig.basis;
            end

            if route.requiresSubsolverSilence
                obj.silenceSubsolvers(solver);
            end
        end

        function silenceSubsolvers(~, solver)
            % Disable reports from nested solvers
            %
            % Args:
            %     solver (handle): Solver object that may own nested solvers
            subsolverNames = {'equilibriumSolver', 'shockSolver', 'jumpConditionsSolver'};

            for i = 1:numel(subsolverNames)
                name = subsolverNames{i};

                if isprop(solver, name) && isobject(solver.(name)) && isprop(solver.(name), 'FLAG_RESULTS')
                    solver.(name).FLAG_RESULTS = false;
                end
            end
        end

        function configureShockTurbulenceSolver(~, solver, solverProblemType)
            % Set properties of the ShockTurbulenceSolver object
            %
            % Args:
            %     solver (ShockTurbulenceSolver): Shock-turbulence solver object
            %     solverProblemType (char): Shock-turbulence model identifier
            solver.setShockTurbulenceModel(solverProblemType);

            % Shock-turbulence is based on planar incident shock jump conditions
            if isprop(solver, 'shockSolver') && isobject(solver.shockSolver) ...
                    && isprop(solver.shockSolver, 'problemType')
                solver.shockSolver.problemType = 'SHOCK_I';
            end

            if isprop(solver, 'jumpConditionsSolver') && isobject(solver.jumpConditionsSolver) ...
                    && isprop(solver.jumpConditionsSolver, 'shockSolver') ...
                    && isobject(solver.jumpConditionsSolver.shockSolver) ...
                    && isprop(solver.jumpConditionsSolver.shockSolver, 'problemType')
                solver.jumpConditionsSolver.shockSolver.problemType = 'SHOCK_I';
            end
        end

        function configureReport(~, solver, options)
            % Configure result report and plotting options
            %
            % Args:
            %     solver (handle): Solver object used for the current problem
            %     options (struct): Report and display options from the setup
            if ~isprop(solver, 'plotConfig') || isempty(solver.plotConfig)
                return
            end

            solver.plotConfig.mintolDisplay = options.displaySpeciesTolerance;
            solver.plotConfig.displaySpecies = options.displaySpecies;
        end

        function solverOutputs = solveWithProblemOptions(obj, solver, route, problem, layout)
            % Solve a problem with scoped solver options from the setup
            %
            % Args:
            %     solver (handle): Solver object used for the current problem
            %     route (struct): Solver route metadata from AppProblemCatalog
            %     problem (struct): Solver-ready problem setup
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %
            % Returns:
            %     solverOutputs (cell): Raw solver outputs arranged for AppResults
            previousModels = obj.applyFrozenChemistryModel(solver, problem);
            cleanup = onCleanup(@() obj.restoreCaloricGasModels(previousModels));
            solverOutputs = obj.solveProblem(solver, route, problem, layout);
        end

        function previousModels = applyFrozenChemistryModel(obj, solver, problem)
            % Apply the thermally perfect model for frozen-chemistry setups
            %
            % Args:
            %     solver (handle): Solver object used for the current problem
            %     problem (struct): Solver-ready problem setup
            %
            % Returns:
            %     previousModels (struct): Solver handles and previous caloric gas models
            previousModels = struct('solver', {}, 'caloricGasModel', {});

            if ~obj.problemUsesFrozenChemistry(problem)
                return
            end

            previousModels = obj.setThermallyPerfectModel(previousModels, solver);

            if isprop(solver, 'equilibriumSolver') && isobject(solver.equilibriumSolver)
                previousModels = obj.setThermallyPerfectModel(previousModels, solver.equilibriumSolver);
            end
        end

        function value = problemUsesFrozenChemistry(~, problem)
            % Check whether a problem setup requires frozen chemistry
            %
            % Args:
            %     problem (struct): Solver-ready problem setup
            %
            % Returns:
            %     value (logical): True when frozen chemistry is active
            value = isfield(problem, 'flags') ...
                && isfield(problem.flags, 'frozenChemistry') ...
                && problem.flags.frozenChemistry;
        end

        function previousModels = setThermallyPerfectModel(~, previousModels, solver)
            % Set one solver to thermally perfect and store its previous model
            %
            % Args:
            %     previousModels (struct): Previously stored caloric gas models
            %     solver (handle): Solver object with a caloricGasModel property
            %
            % Returns:
            %     previousModels (struct): Updated caloric gas model restoration data
            if isempty(solver) || ~isprop(solver, 'caloricGasModel')
                return
            end

            previousModels(end + 1).solver = solver;
            previousModels(end).caloricGasModel = solver.caloricGasModel;
            solver.caloricGasModel = solver.caloricGasModel.setThermallyPerfect();
        end

        function restoreCaloricGasModels(~, previousModels)
            % Restore caloric gas models after a scoped solver option change
            %
            % Args:
            %     previousModels (struct): Solver handles and previous caloric gas models
            for i = numel(previousModels):-1:1
                solver = previousModels(i).solver;

                if isvalid(solver) && isprop(solver, 'caloricGasModel')
                    solver.caloricGasModel = previousModels(i).caloricGasModel;
                end
            end
        end

        function solverOutputs = solveProblem(obj, solver, route, problem, layout)
            % Dispatch a solver call using the configured route
            %
            % Args:
            %     solver (handle): Solver object used for the current problem
            %     route (struct): Solver route metadata from AppProblemCatalog
            %     problem (struct): Solver-ready problem setup
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %
            % Returns:
            %     solverOutputs (cell): Raw solver outputs arranged for AppResults
            switch route.solver
                case 'equilibriumSolver'
                    solverOutputs = obj.solveEquilibriumProblem(solver, problem);
                case 'shockTurbulenceSolver'
                    solverOutputs = cell(1, obj.outputCount(layout, problem.problemType));
                    [solverOutputs{:}] = solver.solve(problem.mixArrayReactants);
                otherwise
                    solverOutputs = cell(1, obj.outputCount(layout, problem.problemType));
                    [solverOutputs{:}] = solver.solveArray(problem.mixArrayReactants);
            end
        end

        function solverOutputs = solveEquilibriumProblem(~, solver, problem)
            % Solve an equilibrium problem from the products mixture array
            %
            % Args:
            %     solver (EquilibriumSolver): Equilibrium solver object
            %     problem (struct): Solver-ready problem setup
            %
            % Returns:
            %     solverOutputs (cell): Reactants and products mixture arrays
            solver.solveArray(problem.mixArrayProducts);
            solverOutputs = {problem.mixArrayReactants, problem.mixArrayProducts};
        end

        function outputs = selectReportOutputs(~, route, solverOutputs)
            % Select the solver outputs used for automatic reports
            %
            % Args:
            %     route (struct): Solver route metadata from AppProblemCatalog
            %     solverOutputs (cell): Raw solver outputs
            %
            % Returns:
            %     outputs (cell): Outputs passed to report generation
            outputs = solverOutputs;

            if strcmp(route.solver, 'equilibriumSolver') && numel(solverOutputs) >= 2
                outputs = solverOutputs(2);
            end
        end

        function results = attachSetupMetadata(~, results, setup)
            % Attach setup metadata to semantic GUI results
            %
            % Args:
            %     results (struct): Semantic result cases from AppResults
            %     setup (struct): Canonical problem setup
            %
            % Returns:
            %     results (struct): Result cases with setup metadata
            for i = 1:numel(results)
                results(i).Reactants = combustiontoolbox.gui.AppSolver.textValue(setup.reactants.selection, 'Custom');
                results(i).Products = combustiontoolbox.gui.AppSolver.productLabel(setup);
                results(i).UITable_R_Data = setup.reactants.tableData;
                results(i).setup = setup;
                results(i).caseIndex = i;

                selectedMixture = results(i).mix2;
                results(i).listSpecies = selectedMixture.chemicalSystem.listSpecies;
                results(i).listProducts = selectedMixture.chemicalSystem.listProducts;
            end
        end

        function count = outputCount(obj, layout, variant)
            % Return the number of outputs expected for a result variant
            %
            % Args:
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %     variant (char): Solver output variant
            %
            % Returns:
            %     count (float): Number of output states
            outputStates = obj.outputStatesForVariant(layout, variant);
            count = numel(outputStates);
        end

        function outputStates = outputStatesForVariant(~, layout, variant)
            % Return output states for the selected result variant
            %
            % Args:
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %     variant (char): Solver output variant
            %
            % Returns:
            %     outputStates (struct): Output state descriptors
            outputStates = layout.outputStates;

            if ~isfield(layout, 'variantOutputStates')
                return
            end

            variant = char(variant);
            variantNames = fieldnames(layout.variantOutputStates);

            for i = 1:numel(variantNames)
                variantName = variantNames{i};

                if strcmpi(variant, variantName) || endsWith(upper(variant), ['_', upper(variantName)])
                    outputStates = layout.variantOutputStates.(variantName);
                    return
                end
            end
        end

        function value = solverProblemType(~, problemType, route)
            % Convert GUI problem type to the solver-level problem type
            %
            % Args:
            %     problemType (char): GUI problem type identifier
            %     route (struct): Solver route metadata from AppProblemCatalog
            %
            % Returns:
            %     value (char): Solver-level problem type identifier
            if strcmp(route.solver, 'shockTurbulenceSolver')
                value = lower(char(route.solverProblemType));
                return
            end

            value = char(problemType);
            value = strrep(value, '_BETA', '');
            value = strrep(value, '_THETA', '');
            value = strrep(value, '_ARATIO', '');
        end
    end

    methods (Static, Access = private)
        function value = textValue(inputValue, defaultValue)
            % Convert a value to text with a fallback
            %
            % Args:
            %     inputValue: Value to convert
            %     defaultValue (char): Fallback text
            %
            % Returns:
            %     value (char): Text value
            if isempty(inputValue)
                value = defaultValue;
                return
            end

            if ischar(inputValue)
                value = inputValue;
            elseif isstring(inputValue)
                value = char(inputValue);
            elseif isnumeric(inputValue) || islogical(inputValue)
                value = num2str(inputValue);
            else
                value = defaultValue;
            end

            if isempty(value)
                value = defaultValue;
            end
        end

        function value = productLabel(setup)
            % Return the product-set label shown in result tree nodes
            %
            % Args:
            %     setup (struct): Canonical problem setup
            %
            % Returns:
            %     value (char): Product-set label
            if isfield(setup, 'flags') && isfield(setup.flags, 'frozenChemistry') ...
                    && setup.flags.frozenChemistry
                value = 'Frozen';
                return
            end

            value = combustiontoolbox.gui.AppSolver.textValue(setup.products.selection, 'All');
        end
    end
end

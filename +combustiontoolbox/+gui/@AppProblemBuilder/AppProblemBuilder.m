classdef AppProblemBuilder < handle
    % Builds canonical and solver-ready problem setups from GUI inputs
    %
    % Attributes:
    %     session (AppSession): Long-lived GUI services
    %     equilibriumProblems (cell): Equilibrium problem type identifiers
    %     constantVolumeProblems (cell): Problem types using volume state input
    %
    % Examples:
    %     * builder = combustiontoolbox.gui.AppProblemBuilder(session);
    %     * setup = builder.buildSetup(input);
    %     * problem = builder.build(input);

    properties (Access = private)
        session
    end

    properties (Access = private, Constant)
        equilibriumProblems = {'TP', 'HP', 'SP', 'TV', 'EV', 'SV'}
        constantVolumeProblems = {'TV', 'EV', 'SV'}
    end

    methods
        function obj = AppProblemBuilder(session)
            % AppProblemBuilder constructor
            %
            % Args:
            %     session (AppSession): Long-lived GUI services
            %
            % Returns:
            %     obj (AppProblemBuilder): Initialized problem builder
            if nargin
                obj.session = session;
            end
        end

        function problem = build(obj, input, varargin)
            % Build a solver-ready problem setup from an input snapshot
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Optional Args:
            %     mixture (Mixture): Current mixture state from the GUI
            %
            % Returns:
            %     problem (struct): Solver-ready problem setup
            currentMixture = [];

            if nargin > 2
                currentMixture = varargin{1};
            end

            setup = obj.buildSetup(input, currentMixture);
            obj.validateReactantsForBuild(setup.reactants);
            problem = obj.buildBaseProblem(setup.productSpecies);
            problem.setup = setup;
            problem.input = input;
            problem.requestedProblemType = setup.requestedProblemType;
            problem.problemType = setup.problemType;
            problem.productSpecies = setup.productSpecies;
            problem.flags = setup.flags;
            problem.additionalInputsReactants = setup.additionalInputsReactants;
            problem.additionalInputsProducts = setup.additionalInputsProducts;
            problem.reactantStateProperties = setup.reactantStateProperties;
            problem.productStateProperties = setup.productStateProperties;

            % Build the reactant mixture and the first CT mixture array
            problem.mixture = obj.buildMixture(problem.chemicalSystem, input, currentMixture);
            reactantStateInputs = obj.statePropertiesToInputs(problem.reactantStateProperties);
            problem.mixArrayReactants = setProperties(problem.mixture, ...
                reactantStateInputs{:}, problem.additionalInputsReactants{:});
            obj.applySpeciesTemperatures(problem.mixArrayReactants, input.reactantsTable);
            problem = obj.setupProductStateProperties(problem, input);
            problem.setup.productStateProperties = problem.productStateProperties;
            problem.setup.targetState.properties = problem.productStateProperties;

            if obj.isEquilibriumProblem(problem.requestedProblemType)
                productStateInputs = obj.statePropertiesToInputs(problem.productStateProperties);
                problem.mixArrayProducts = setProperties(problem.mixture, ...
                    productStateInputs{:}, problem.additionalInputsProducts{:});
            end
        end

        function setup = buildSetup(obj, input, varargin)
            % Build a canonical problem setup from GUI input data
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Optional Args:
            %     mixture (Mixture): Current mixture state from the GUI
            %
            % Returns:
            %     setup (struct): Solver-independent setup description
            currentMixture = [];

            if nargin > 2
                currentMixture = varargin{1};
            end

            reactants = obj.reactantsSetup(input, currentMixture);
            productSpecies = obj.resolveProductSpecies(input);
            setup = struct();
            setup.requestedProblemType = char(input.problemType);
            setup.problemType = char(input.problemType);
            setup.constraint = obj.constraintName(input.problemType);
            setup.family = obj.problemFamily(input.problemType);
            setup.reactants = reactants;
            setup.products = struct( ...
                'selection', input.products, ...
                'type', obj.productSelectionType(input, reactants, productSpecies), ...
                'species', {productSpecies});
            setup.productSpecies = productSpecies;
            setup.flags = obj.defaultFlags(input, reactants);
            setup.options = obj.optionsSetup(input);
            setup.additionalInputsReactants = obj.equivalenceRatioInputs(input, setup.reactants);
            setup.additionalInputsProducts = obj.equivalenceRatioInputs(input, setup.reactants);
            setup.additionalProperties = input.additionalProperties;
            setup.expressions = input.exportExpressions();
            setup.reactantStateProperties = obj.makeStateProperties( ...
                input, input.temperatureReactants, input.pressureReactants);
            setup.productStateProperties = setup.reactantStateProperties;
            setup.initialState = struct('properties', setup.reactantStateProperties);
            setup.targetState = struct('properties', setup.productStateProperties);
            setup = obj.applyProblemTypeInputs(setup, input);
        end

        function script = exportScript(~, setup)
            % Export a canonical problem setup to a MATLAB script
            %
            % Args:
            %     setup (struct): Canonical problem setup
            %
            % Returns:
            %     script (char): Runnable MATLAB script text
            exporter = combustiontoolbox.gui.AppProblemScriptExporter();
            script = exporter.export(setup);
        end

    end

    methods (Access = private)
        function setup = reactantsSetup(obj, input, varargin)
            % Build the canonical reactants setup from the GUI table
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Optional Args:
            %     mixture (Mixture): Current mixture state from the GUI
            %
            % Returns:
            %     setup (struct): Reactant species, amounts, roles, and temperatures
            setup = struct( ...
                'selection', input.reactants, ...
                'species', {{}}, ...
                'amount', [], ...
                'role', {{}}, ...
                'temperature', [], ...
                'tableData', {{}});
            currentMixture = [];

            if nargin > 2
                currentMixture = varargin{1};
            end

            data = input.reactantsTable;

            if isempty(data)
                return
            end

            setup.species = data(:, 1)';
            setup.amount = obj.cellColumnToRowVector(data(:, 2));
            setup.role = data(:, 4)';
            setup.tableData = data;

            if size(data, 2) >= 5
                setup.temperature = obj.cellColumnToRowVector(data(:, 5));
            end

            if ~isempty(input.equivalenceRatio) && obj.hasMixtureComposition(currentMixture)
                setup.amount = obj.referenceAmountsForSpecies( ...
                    setup.species, setup.role, currentMixture, setup.amount);
            end
        end

        function validateReactantsForBuild(~, reactants)
            % Validate that a problem can be built from the reactants setup
            %
            % Args:
            %     reactants (struct): Canonical reactants setup
            if isempty(reactants.species)
                error('AppProblemBuilder:MissingReactants', ...
                    'Select at least one reactant before solving the problem.');
            end
        end

        function options = optionsSetup(~, input)
            % Build report and display options from GUI inputs
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     options (struct): Report and display options
            printResults = false;

            if isfield(input.flags, 'printResults')
                printResults = input.flags.printResults;
            end

            options = struct( ...
                'printResults', printResults, ...
                'reportType', input.reportType, ...
                'displaySpecies', {input.displaySpecies}, ...
                'displaySpeciesTolerance', input.displaySpeciesTolerance);
        end

        function value = resolveProductSpecies(~, input)
            % Resolve the product species representation used by ChemicalSystem
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     value (cell | char): Product species or predefined product-set keyword
            value = input.productSpecies;

            if strcmpi(input.products, 'complete reaction')
                value = 'complete';
            end
        end

        function value = productSelectionType(obj, input, reactants, productSpecies)
            % Classify the selected product species source
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     reactants (struct): Canonical reactants setup
            %     productSpecies (cell | char): Product species resolved for the setup
            %
            % Returns:
            %     value (char): Product selection type
            if strcmpi(input.products, 'complete reaction')
                value = 'complete';
            elseif obj.hasOnlyInertReactants(reactants)
                value = 'frozen';
            elseif isempty(input.products) ...
                    && obj.matchesAutomaticProductSpecies(input, reactants, productSpecies)
                value = 'auto';
            else
                value = 'explicit';
            end
        end

        function value = matchesAutomaticProductSpecies(obj, input, reactants, productSpecies)
            % Check whether the current product list matches the automatic list
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     reactants (struct): Canonical reactants setup
            %     productSpecies (cell | char): Product species resolved for the setup
            %
            % Returns:
            %     value (logical): True when the product species list is automatic
            value = false;

            if isempty(productSpecies) || ischar(productSpecies) ...
                    || isempty(obj.session) || isempty(obj.session.database)
                return
            end

            automaticSpecies = obj.automaticProductSpecies(input, reactants);
            value = isequal(productSpecies(:)', automaticSpecies(:)');
        end

        function value = automaticProductSpecies(obj, input, reactants)
            % Compute the automatic product species for the reactant set
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     reactants (struct): Canonical reactants setup
            %
            % Returns:
            %     value (cell): Automatically detected product species
            value = {};

            if isempty(reactants.species)
                return
            end

            reactiveSpecies = obj.reactiveProductSourceSpecies(reactants);

            if isempty(reactiveSpecies)
                return
            end

            chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.session.database);
            value = findProducts(chemicalSystem, reactiveSpecies, ...
                'flag_ion', obj.inputFlag(input, 'ionizedSpecies', false));
        end

        function value = reactiveProductSourceSpecies(~, reactants)
            % Return reactant species that define automatic equilibrium products
            %
            % Args:
            %     reactants (struct): Canonical reactants setup
            %
            % Returns:
            %     value (cell): Fuel and oxidizer species
            role = reactants.role;
            index = strcmpi(role, 'Fuel') | strcmpi(role, 'Oxidizer');
            value = reactants.species(index);
        end

        function value = hasOnlyInertReactants(~, reactants)
            % Check whether all reactants are inert species
            %
            % Args:
            %     reactants (struct): Canonical reactants setup
            %
            % Returns:
            %     value (logical): True when every reactant role is inert
            value = ~isempty(reactants.species) ...
                && all(strcmpi(reactants.role, 'Inert'));
        end

        function flags = defaultFlags(obj, input, reactants)
            % Build solver and setup flags from GUI inputs
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     reactants (struct): Canonical reactants setup
            %
            % Returns:
            %     flags (struct): Solver and setup flags
            printResults = obj.inputFlag(input, 'printResults', false);
            frozenChemistry = obj.inputFlag(input, 'frozenChemistry', false) ...
                || obj.hasOnlyInertReactants(reactants);

            flags = struct( ...
                'hasEquivalenceRatio', ~isempty(input.equivalenceRatio), ...
                'frozenChemistry', frozenChemistry, ...
                'ionizedSpecies', obj.inputFlag(input, 'ionizedSpecies', false), ...
                'idealAir', obj.inputFlag(input, 'idealAir', false), ...
                'infiniteAreaChamber', obj.inputFlag(input, 'infiniteAreaChamber', true), ...
                'printResults', printResults, ...
                'beta', false, ...
                'theta', false, ...
                'areaRatio', false, ...
                'rocketSubsonic', []);
        end

        function value = inputFlag(~, input, name, defaultValue)
            % Read a boolean flag from an AppInput snapshot
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     name (char): Flag name
            %     defaultValue: Fallback value
            %
            % Returns:
            %     value: Flag value
            value = defaultValue;

            if isfield(input.flags, name)
                value = input.flags.(name);
            end
        end

        function inputs = equivalenceRatioInputs(~, input, reactants) %#ok<INUSD>
            % Build name-value inputs for equivalence-ratio dependent calls
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     reactants (struct): Canonical reactants setup
            %
            % Returns:
            %     inputs (cell): Name-value inputs for setProperties
            inputs = {};

            if ~isempty(input.equivalenceRatio)
                inputs = {'equivalenceRatio', input.equivalenceRatio};
            end
        end

        function problem = buildBaseProblem(obj, listSpecies)
            % Build the base problem structure and chemical system
            %
            % Optional Args:
            %     listSpecies (cell | char): Product species or predefined product-set keyword
            %
            % Returns:
            %     problem (struct): Base problem structure
            if nargin < 2
                listSpecies = {};
            end

            problem = struct();
            problem.chemicalSystem = obj.buildChemicalSystem(listSpecies);
            problem.mixture = [];
        end

        function properties = makeStateProperties(obj, input, temperature, mechanicalValue)
            % Build state-property descriptors for a thermodynamic state
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     temperature (float): Temperature input [K]
            %     mechanicalValue (float): Pressure [bar] or volume [m3/kg]
            %
            % Returns:
            %     properties (struct): State-property descriptors
            mechanicalName = 'pressure';

            if obj.isConstantVolumeProblem(input.problemType)
                mechanicalName = 'volume';
            end

            properties = [
                struct('name', 'temperature', 'value', temperature), ...
                struct('name', mechanicalName, 'value', mechanicalValue)];
        end

        function value = isConstantVolumeProblem(obj, problemType)
            % Check whether a problem type uses volume as mechanical state
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (logical): True for constant-volume problem types
            value = any(strcmpi(char(problemType), obj.constantVolumeProblems));
        end

        function value = isEquilibriumProblem(obj, problemType)
            % Check whether a problem type belongs to equilibrium problems
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (logical): True for equilibrium problem types
            value = any(strcmpi(char(problemType), obj.equilibriumProblems));
        end

        function problem = applyProblemTypeInputs(obj, problem, input)
            % Apply inputs and solver variants that depend on the problem type
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Setup with problem-specific inputs and variants
            problem.problemType = char(input.problemType);

            switch problem.problemType
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR'}
                    problem = obj.setupShockProblem(problem, input);
                case {'SHOCK_OBLIQUE', 'SHOCK_OBLIQUE_R', 'SHOCK_POLAR_R'}
                    problem = obj.setupShockObliqueProblem(problem, input);
                case {'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', 'DET_OVERDRIVEN_R', ...
                        'DET_UNDERDRIVEN_R', 'DET_POLAR'}
                    problem = obj.setupDetonationProblem(problem, input);
                case {'DET_OBLIQUE', 'DET_POLAR_R'}
                    problem = obj.setupDetonationObliqueProblem(problem, input);
                case 'ROCKET'
                    problem = obj.setupRocketProblem(problem, input);
                case {'SHOCKTURBULENCE_VORTICAL', 'SHOCKTURBULENCE_ACOUSTIC'}
                    problem = obj.setupShockTurbulenceProblem(problem, input);
                case 'SHOCKTURBULENCE_VORTICAL_ENTROPIC'
                    problem = obj.setupShockTurbulenceVorticalEntropicProblem(problem, input);
                case 'SHOCKTURBULENCE_COMPRESSIBLE'
                    problem = obj.setupShockTurbulenceCompressibleProblem(problem, input);
            end
        end

        function problem = setupShockProblem(obj, problem, input)
            % Add incident shock inputs to the problem setup
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Shock problem setup
            problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                {'mach', obj.additionalPropertyValue(input, 'PR4')}];
        end

        function problem = setupShockObliqueProblem(obj, problem, input)
            % Select beta or theta from the active angle field
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Oblique shock problem setup
            mach = obj.additionalPropertyValue(input, 'PR4');
            beta = obj.additionalPropertyValue(input, 'PR5');
            theta = obj.additionalPropertyValue(input, 'PP5');
            problem = obj.applyAngleVariant(problem, {'mach', mach}, beta, theta);
        end

        function problem = setupDetonationProblem(obj, problem, input)
            % Add detonation drive-factor inputs to the problem setup
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Detonation problem setup
            problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                {'driveFactor', obj.additionalPropertyValue(input, 'PR3')}];
        end

        function problem = setupDetonationObliqueProblem(obj, problem, input)
            % Select beta or theta for oblique detonation variants
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Oblique detonation problem setup
            driveFactor = obj.additionalPropertyValue(input, 'PR3');
            beta = obj.additionalPropertyValue(input, 'PR4');
            theta = obj.additionalPropertyValue(input, 'PP4');
            problem = obj.applyAngleVariant(problem, {'driveFactor', driveFactor}, beta, theta);
        end

        function problem = setupRocketProblem(obj, problem, input)
            % Chamber model is encoded in the solver problem type
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Rocket problem setup
            infiniteAreaChamber = true;

            if isfield(input.flags, 'infiniteAreaChamber')
                infiniteAreaChamber = input.flags.infiniteAreaChamber;
            end

            if infiniteAreaChamber
                problem.problemType = 'ROCKET_IAC';
            else
                problem.problemType = 'ROCKET_FAC';
                areaRatioChamber = input.temperatureProducts;
                assert(~isempty(areaRatioChamber), 'Specify the area chamber / area throat ratio.');
                problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                    {'areaRatioChamber', areaRatioChamber}];
            end

            % Optional area ratio selects the nozzle branch
            subsonicAreaRatio = obj.additionalPropertyValue(input, 'PR3');
            supersonicAreaRatio = obj.additionalPropertyValue(input, 'PP3');

            if ~isempty(subsonicAreaRatio)
                problem.flags.areaRatio = true;
                problem.flags.rocketSubsonic = true;
                problem.problemType = [problem.problemType, '_ARATIO'];
                problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                    {'areaRatio', subsonicAreaRatio}];
            elseif ~isempty(supersonicAreaRatio)
                problem.flags.areaRatio = true;
                problem.flags.rocketSubsonic = false;
                problem.problemType = [problem.problemType, '_ARATIO'];
                problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                    {'areaRatio', supersonicAreaRatio}];
            end
        end

        function problem = setupShockTurbulenceProblem(obj, problem, input)
            % Add base shock-turbulence inputs to the problem setup
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Shock-turbulence problem setup
            problem.productStateProperties(1).value = 0;
            problem.productStateProperties(2).value = 0;
            problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                {'mach', obj.additionalPropertyValue(input, 'PR4')}];
        end

        function problem = setupShockTurbulenceVorticalEntropicProblem(obj, problem, input)
            % Add vortical-entropic shock-turbulence inputs to the problem setup
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Vortical-entropic shock-turbulence setup
            problem.productStateProperties(1).value = 0;
            problem.productStateProperties(2).value = input.temperatureProducts;
            problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                {'mach', obj.additionalPropertyValue(input, 'PR4'), ...
                'chi', input.temperatureProducts}];
        end

        function problem = setupShockTurbulenceCompressibleProblem(obj, problem, input)
            % Add compressible shock-turbulence inputs to the problem setup
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Compressible shock-turbulence setup
            problem.productStateProperties(1).value = input.temperatureProducts;
            problem.productStateProperties(2).value = input.pressureProducts;
            problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                {'mach', obj.additionalPropertyValue(input, 'PR4'), ...
                'eta', input.temperatureProducts, ...
                'chi', input.pressureProducts, ...
                'etaVorticity', obj.additionalPropertyValue(input, 'PP6')}];
        end

        function problem = applyAngleVariant(~, problem, leadingInputs, beta, theta)
            % Append beta or theta inputs and tag the solver variant
            %
            % Args:
            %     problem (struct): Canonical problem setup
            %     leadingInputs (cell): Name-value inputs required before the angle
            %     beta (float): Wave angle input [deg]
            %     theta (float): Deflection angle input [deg]
            %
            % Returns:
            %     problem (struct): Problem setup with angle variant metadata
            if ~isempty(beta)
                problem.flags.beta = true;
                problem.problemType = [problem.problemType, '_BETA'];
                problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                    leadingInputs, {'beta', beta}];
            else
                problem.flags.theta = true;
                problem.problemType = [problem.problemType, '_THETA'];
                problem.additionalInputsReactants = [problem.additionalInputsReactants, ...
                    leadingInputs, {'theta', theta}];
            end
        end

        function inputs = statePropertiesToInputs(~, stateProperties)
            % Convert state property structs to CT name-value inputs
            %
            % Args:
            %     stateProperties (struct): State-property descriptors
            %
            % Returns:
            %     inputs (cell): Name-value inputs for setProperties
            inputs = cell(1, 2 * numel(stateProperties));

            for i = 1:numel(stateProperties)
                inputs{2 * i - 1} = stateProperties(i).name;
                inputs{2 * i} = stateProperties(i).value;
            end
        end

        function value = additionalPropertyValue(~, input, name)
            % Read one additional GUI property by name
            %
            % Args:
            %     input (AppInput): Typed GUI input snapshot
            %     name (char): Additional property field name
            %
            % Returns:
            %     value: Additional property value or empty
            value = [];

            if isfield(input.additionalProperties, name)
                value = input.additionalProperties.(name);
            end
        end

        function mixture = buildMixture(obj, chemicalSystem, input, currentMixture)
            % Build a mixture using current composition or table data
            %
            % Args:
            %     chemicalSystem (ChemicalSystem): Chemical system for the mixture
            %     input (AppInput): Typed GUI input snapshot
            %     currentMixture (Mixture): Current GUI mixture state
            %
            % Returns:
            %     mixture (Mixture): Mixture with reactant composition
            mixture = combustiontoolbox.core.Mixture(chemicalSystem);

            if ~isempty(input.equivalenceRatio) && obj.hasMixtureComposition(currentMixture)
                obj.copyMixtureComposition(mixture, currentMixture);
            elseif ~isempty(input.reactantsTable)
                obj.applyReactantsTableComposition(mixture, input.reactantsTable);
            elseif obj.hasMixtureComposition(currentMixture)
                obj.copyMixtureComposition(mixture, currentMixture);
            end
        end

        function value = hasMixtureComposition(~, mixture)
            % Check whether a mixture contains any composition data
            %
            % Args:
            %     mixture (Mixture): Mixture object to inspect
            %
            % Returns:
            %     value (logical): True when fuel, oxidizer, or inert species exist
            value = false;

            if ~isobject(mixture)
                return
            end

            speciesProperties = {'listSpeciesFuel', 'listSpeciesOxidizer', 'listSpeciesInert'};

            for i = 1:numel(speciesProperties)
                name = speciesProperties{i};

                if isprop(mixture, name) && ~isempty(mixture.(name))
                    value = true;
                    return
                end
            end
        end

        function problem = setupProductStateProperties(obj, problem, input)
            % Build product state-property descriptors for equilibrium problems
            %
            % Args:
            %     problem (struct): Solver-ready problem setup
            %     input (AppInput): Typed GUI input snapshot
            %
            % Returns:
            %     problem (struct): Problem setup with product state properties
            productPressure = input.pressureProducts;

            if isempty(productPressure)
                productPressure = input.pressureReactants;
            end

            switch char(input.problemType)
                case {'TP', 'TV'}
                    problem.productStateProperties = obj.makeStateProperties( ...
                        input, input.temperatureProducts, input.pressureProducts);
                case 'HP'
                    problem.productStateProperties = obj.makeStateProperties( ...
                        input, input.temperatureReactants, productPressure);
                case 'SP'
                    entropy = [problem.mixArrayReactants.sSpecific];
                    problem.productStateProperties = [
                        struct('name', 'entropySpecific', 'value', entropy), ...
                        struct('name', 'pressure', 'value', productPressure)];
                case 'SV'
                    entropy = [problem.mixArrayReactants.sSpecific];
                    volumeRatio = obj.additionalPropertyValue(input, 'PP4');

                    if isempty(volumeRatio)
                        volumeRatio = 1;
                    end

                    productVolume = input.pressureReactants .* volumeRatio;
                    problem.productStateProperties = [
                        struct('name', 'entropySpecific', 'value', entropy), ...
                        struct('name', 'volume', 'value', productVolume)];
                otherwise
                    problem.productStateProperties = problem.reactantStateProperties;
            end
        end

        function value = constraintName(~, problemType)
            % Return the thermodynamic constraint name for equilibrium problems
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (char): Constraint name
            switch upper(char(problemType))
                case 'TP'
                    value = 'temperaturePressure';
                case 'TV'
                    value = 'temperatureVolume';
                case 'HP'
                    value = 'enthalpyPressure';
                case 'EV'
                    value = 'internalEnergyVolume';
                case 'SP'
                    value = 'entropyPressure';
                case 'SV'
                    value = 'entropyVolume';
                otherwise
                    value = '';
            end
        end

        function value = problemFamily(~, problemType)
            % Return the problem family for a problem type
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (char): Problem family name
            problemType = upper(char(problemType));

            if any(strcmp(problemType, {'TP', 'HP', 'SP', 'TV', 'EV', 'SV'}))
                value = 'equilibrium';
            elseif startsWith(problemType, 'SHOCKTURBULENCE')
                value = 'shockTurbulence';
            elseif startsWith(problemType, 'SHOCK')
                value = 'shock';
            elseif startsWith(problemType, 'DET')
                value = 'detonation';
            elseif startsWith(problemType, 'ROCKET')
                value = 'rocket';
            else
                value = 'unknown';
            end
        end

        function copyMixtureComposition(~, mixture, source)
            % Copy fuel, oxidizer, and inert composition from an existing mixture
            %
            % Args:
            %     mixture (Mixture): Destination mixture
            %     source (Mixture): Source mixture
            if ~isempty(source.listSpeciesFuel)
                set(mixture, source.listSpeciesFuel, 'fuel', source.molesFuel);
            end

            if ~isempty(source.listSpeciesOxidizer)
                ratioOxidizer = source.ratioOxidizer;

                if isempty(ratioOxidizer)
                    ratioOxidizer = source.molesOxidizer;
                end

                set(mixture, source.listSpeciesOxidizer, 'oxidizer', ratioOxidizer);
                mixture.ratioOxidizer = ratioOxidizer;
            end

            if ~isempty(source.listSpeciesInert)
                set(mixture, source.listSpeciesInert, 'inert', source.molesInert);
            end
        end

        function amounts = referenceAmountsForSpecies(~, species, roles, mixture, fallbackAmounts)
            % Return composition amounts from the current mixture when possible
            %
            % Args:
            %     species (cell): Species names from the reactants table
            %     roles (cell): Species roles from the reactants table
            %     mixture (Mixture): Current GUI mixture state
            %     fallbackAmounts (float): Amounts from the reactants table
            %
            % Returns:
            %     amounts (float): Reference amounts aligned with species
            amounts = fallbackAmounts;
            oxidizerAmounts = mixture.ratioOxidizer;

            if isempty(oxidizerAmounts)
                oxidizerAmounts = mixture.molesOxidizer;
            end

            groups = { ...
                'Fuel', mixture.listSpeciesFuel, mixture.molesFuel; ...
                'Oxidizer', mixture.listSpeciesOxidizer, oxidizerAmounts; ...
                'Inert', mixture.listSpeciesInert, mixture.molesInert};

            for i = 1:size(groups, 1)
                roleName = groups{i, 1};
                groupSpecies = groups{i, 2};
                groupAmounts = groups{i, 3};
                index = strcmpi(roles, roleName);

                if ~any(index) || isempty(groupSpecies)
                    continue
                end

                [found, sourceIndex] = ismember(species(index), groupSpecies);
                targetIndex = find(index);
                targetIndex = targetIndex(found);
                amounts(targetIndex) = groupAmounts(sourceIndex(found));
            end
        end

        function applyReactantsTableComposition(obj, mixture, data)
            % Build reactant composition from table data
            %
            % Args:
            %     mixture (Mixture): Destination mixture
            %     data (cell): Reactants table data
            if isempty(data)
                return
            end

            listSpecies = data(:, 1)';
            moles = obj.cellColumnToRowVector(data(:, 2));
            typeSpecies = data(:, 4);
            indexFuel = contains(typeSpecies, 'Fuel');
            indexOxidizer = contains(typeSpecies, 'Oxidizer');
            indexInert = contains(typeSpecies, 'Inert');

            if any(indexFuel)
                set(mixture, listSpecies(indexFuel), 'fuel', moles(indexFuel));
            end

            if any(indexOxidizer)
                set(mixture, listSpecies(indexOxidizer), 'oxidizer', moles(indexOxidizer));
                mixture.ratioOxidizer = moles(indexOxidizer);
            end

            if any(indexInert)
                set(mixture, listSpecies(indexInert), 'inert', moles(indexInert));
            end
        end

        function applySpeciesTemperatures(obj, mixtures, data)
            % Apply species temperatures from the reactants table
            %
            % Args:
            %     mixtures (Mixture): Mixture array to update
            %     data (cell): Reactants table data
            if isempty(data) || size(data, 2) < 5
                return
            end

            speciesTemperatures = obj.cellColumnToRowVector(data(:, 5));

            if obj.isBulkTemperatureSweep(mixtures) && obj.isUniformValue(speciesTemperatures)
                return
            end

            for i = 1:numel(mixtures)
                mixture = mixtures(i);
                [~, index] = ismember(mixture.listSpecies, data(:, 1));

                if isempty(index) || any(index == 0)
                    continue
                end

                temperatures = speciesTemperatures(index);

                if all(abs(temperatures - mixture.T) < 1e-12)
                    continue
                end

                setTemperatureSpecies(mixture, temperatures);
                updateThermodynamics(mixture);
            end
        end

        function value = isBulkTemperatureSweep(~, mixtures)
            % Check whether a mixture array was generated from bulk temperature
            %
            % Args:
            %     mixtures (Mixture): Mixture array to inspect
            %
            % Returns:
            %     value (logical): True for parametric bulk temperature studies
            value = numel(mixtures) > 1 ...
                && isprop(mixtures(1), 'rangeName') ...
                && strcmpi(mixtures(1).rangeName, 'T');
        end

        function value = isUniformValue(~, values)
            % Check whether all numeric values are equal within tolerance
            %
            % Args:
            %     values (float): Values to inspect
            %
            % Returns:
            %     value (logical): True when values are empty or uniform
            value = isempty(values) || all(abs(values - values(1)) < 1e-12);
        end

        function chemicalSystem = buildChemicalSystem(obj, listSpecies)
            % Build the core chemical system for a species list
            %
            % Args:
            %     listSpecies (cell | char): Species list or catalog keyword
            %
            % Returns:
            %     chemicalSystem (ChemicalSystem): Core chemical system
            obj.validateSession();

            if isempty(listSpecies)
                chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.session.database);
            else
                chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.session.database, listSpecies);
            end
        end

        function validateSession(obj)
            % Validate services required to create chemical systems
            %
            % Args:
            %     obj (AppProblemBuilder): Problem builder object
            if isempty(obj.session) || isempty(obj.session.database)
                error('AppProblemBuilder:MissingSession', ...
                    'AppProblemBuilder requires an AppSession with a database.');
            end
        end

        function value = cellColumnToRowVector(~, values)
            % Convert table column data to a row vector
            %
            % Args:
            %     values (cell | float): Table column values
            %
            % Returns:
            %     value (float): Row vector
            if iscell(values)
                value = cell2mat(values)';
            else
                value = values';
            end
        end
    end
end

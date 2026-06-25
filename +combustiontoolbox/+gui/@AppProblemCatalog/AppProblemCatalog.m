classdef AppProblemCatalog < handle
    % Defines GUI problem metadata, solver routing, plots, and result layouts.
    %
    % Attributes:
    %     definitions (struct): Problem metadata keyed by valid MATLAB names
    %     order (cell): Stable problem order used by the GUI dropdown
    %
    % Examples:
    %     * catalog = combustiontoolbox.gui.AppProblemCatalog();
    %     * catalog.has('TP')
    %     * [items, itemsData] = catalog.dropdownData();
    %     * definition = catalog.get('ROCKET');

    properties (Access = private)
        definitions
        order
    end

    methods
        function obj = AppProblemCatalog()
            % AppProblemCatalog constructor
            %
            % Returns:
            %     obj (AppProblemCatalog): Initialized problem catalog
            [obj.definitions, obj.order] = obj.createDefinitions();
        end

        function value = list(obj)
            % List all registered problem identifiers in dropdown order
            %
            % Returns:
            %     value (cell): Registered problem identifiers
            value = obj.ids();
        end

        function value = ids(obj)
            % List all registered problem identifiers in dropdown order
            %
            % Returns:
            %     value (cell): Registered problem identifiers
            value = obj.order;
        end

        function value = itemsData(obj)
            % Return dropdown ItemsData values
            %
            % Returns:
            %     value (cell): Problem identifiers used as dropdown values
            value = obj.ids();
        end

        function value = items(obj)
            % Return dropdown Items labels
            %
            % Returns:
            %     value (cell): Human-readable dropdown labels
            value = cell(size(obj.order));

            for i = 1:numel(obj.order)
                definition = obj.get(obj.order{i});
                value{i} = definition.label;
            end

        end

        function [items, itemsData] = dropdownData(obj)
            % Return labels and identifiers for the problem dropdown
            %
            % Returns:
            %     items (cell): Human-readable dropdown labels
            %     itemsData (cell): Problem identifiers used as dropdown values
            items = obj.items();
            itemsData = obj.itemsData();
        end

        function value = all(obj)
            % Return all definitions in dropdown order
            %
            % Returns:
            %     value (cell): Problem definition structs
            value = cell(size(obj.order));

            for i = 1:numel(obj.order)
                value{i} = obj.get(obj.order{i});
            end

        end

        function value = has(obj, problemType)
            % Check whether a problem type is registered
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (logical): True if problemType is registered
            key = obj.normalize(problemType);
            value = isfield(obj.definitions, key);
        end

        function definition = get(obj, problemType)
            % Get a problem definition
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     definition (struct): Problem metadata
            key = obj.normalize(problemType);

            if ~isfield(obj.definitions, key)
                error('AppProblemCatalog:UnknownProblem', ...
                    'Problem type "%s" is not registered in AppProblemCatalog.', problemType);
            end

            definition = obj.definitions.(key);
        end

        function value = byFamily(obj, family)
            % List registered problem identifiers that belong to a family
            %
            % Args:
            %     family (char): Problem family name
            %
            % Returns:
            %     value (cell): Problem identifiers in dropdown order
            value = {};

            for i = 1:numel(obj.order)
                definition = obj.get(obj.order{i});

                if strcmpi(definition.family, family)
                    value{end + 1} = definition.id; %#ok<AGROW>
                end

            end

        end

        function value = families(obj)
            % List registered problem families in first-seen order
            %
            % Returns:
            %     value (cell): Registered family names
            value = {};

            for i = 1:numel(obj.order)
                definition = obj.get(obj.order{i});

                if ~any(strcmp(value, definition.family))
                    value{end + 1} = definition.family; %#ok<AGROW>
                end

            end

        end

        function value = defaults(obj, problemType)
            % Return default GUI input values for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Default values keyed by GUI component name
            definition = obj.get(problemType);
            value = definition.defaultInputs;
        end

        function value = labels(obj, problemType)
            % Return GUI text labels for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Label text keyed by GUI component name
            definition = obj.get(problemType);
            value = definition.labels;
        end

        function value = visibility(obj, problemType)
            % Return visibility metadata for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Logical visibility flags
            definition = obj.get(problemType);
            value = definition.visibility;
        end

        function value = solverRoute(obj, problemType)
            % Return solver routing metadata for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Solver dispatch metadata
            definition = obj.get(problemType);
            value = definition.solverRoute;
        end

        function [properties, basis] = plotProperties(obj, problemType)
            % Return default plot properties for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     properties (cell): Plot property names
            %     basis (cell): Plot property bases
            definition = obj.get(problemType);
            properties = definition.plotConfig.properties;
            basis = definition.plotConfig.basis;
        end

        function value = resultLayout(obj, problemType)
            % Return result layout metadata for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Result state names and layout mode
            definition = obj.get(problemType);
            value = definition.resultLayout;
        end

        function value = status(obj, problemType)
            % Return implementation status for a problem
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (char): Implementation status
            definition = obj.get(problemType);
            value = definition.status;
        end
    end

    methods (Access = private)
        function [definitions, order] = createDefinitions(obj)
            order = { ...
                'TP', ...
                'HP', ...
                'SP', ...
                'TV', ...
                'EV', ...
                'SV', ...
                'SHOCK_I', ...
                'SHOCK_R', ...
                'SHOCK_OBLIQUE', ...
                'SHOCK_OBLIQUE_R', ...
                'SHOCK_POLAR', ...
                'SHOCK_POLAR_R', ...
                'DET', ...
                'DET_OVERDRIVEN', ...
                'DET_UNDERDRIVEN', ...
                'DET_R', ...
                'DET_OVERDRIVEN_R', ...
                'DET_UNDERDRIVEN_R', ...
                'DET_OBLIQUE', ...
                'DET_POLAR', ...
                'ROCKET', ...
                'SHOCKTURBULENCE_VORTICAL', ...
                'SHOCKTURBULENCE_VORTICAL_ENTROPIC', ...
                'SHOCKTURBULENCE_ACOUSTIC', ...
                'SHOCKTURBULENCE_COMPRESSIBLE'};

            dropdownLabels = { ...
                'TP:  Equilibrium composition at defined T and P', ...
                'HP: Adiabatic T and composition at constant P', ...
                'SP:  Isentropic compression/expansion to a specified P', ...
                'TV: Equilibrium composition at defined T and constant V', ...
                'EV: Adiabatic T and composition at constant V', ...
                'SV: Isentropic compresion/expansion to a specified V', ...
                'SHOCK_I: Planar incident shock wave', ...
                'SHOCK_R: Planar reflected shock wave', ...
                'SHOCK_OBLIQUE: Oblique incident shock wave', ...
                'SHOCK_OBLIQUE_R: Oblique reflected shock wave', ...
                'SHOCK_POLAR: Polar shocks', ...
                'SHOCK_POLAR_R: Polar shocks (regular reflections)', ...
                'DET: Chapman-Jouguet Detonation', ...
                'DET_OVERDRIVEN: Overdriven detonation', ...
                'DET_UNDERDRIVEN: Underdriven detonation', ...
                'DET_R: Reflected Chapman-Jouguet detonation', ...
                'DET_OVERDRIVEN_R: Overdriven reflected detonation', ...
                'DET_UNDERDRIVEN_R: Underdriven reflected detonation', ...
                'DET_OBLIQUE: Oblique incident detonation', ...
                'DET_POLAR: Polar detonations', ...
                'ROCKET: Rocket propellant performance', ...
                'SHOCKTURBULENCE_VORTICAL: Vortical shock-turbulence interaction', ...
                'SHOCKTURBULENCE_VORTICAL_ENTROPIC: Vortical-entropic shock-turbulence interaction', ...
                'SHOCKTURBULENCE_ACOUSTIC: Acoustic shock-turbulence interaction', ...
                'SHOCKTURBULENCE_COMPRESSIBLE: Compressible shock-turbulence interaction'};

            families = { ...
                'equilibrium', ...
                'equilibrium', ...
                'equilibrium', ...
                'equilibrium', ...
                'equilibrium', ...
                'equilibrium', ...
                'shock', ...
                'shock', ...
                'shock', ...
                'shock', ...
                'shock', ...
                'shock', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'detonation', ...
                'rocket', ...
                'shockturbulence', ...
                'shockturbulence', ...
                'shockturbulence', ...
                'shockturbulence'};

            definitions = struct();

            for i = 1:numel(order)
                definition = obj.makeDefinition(order{i}, dropdownLabels{i}, families{i});
                definitions = obj.add(definitions, definition);
            end

        end

        function definitions = add(obj, definitions, definition)
            key = obj.normalize(definition.id);
            definitions.(key) = definition;
        end

        function definition = makeDefinition(obj, id, label, family)
            definition = struct( ...
                'id', id, ...
                'label', label, ...
                'family', family, ...
                'status', obj.makeStatus(id), ...
                'defaultInputs', obj.makeDefaults(id), ...
                'labels', obj.makeLabels(id), ...
                'visibility', obj.makeVisibility(id), ...
                'inputRoles', obj.makeInputRoles(id), ...
                'solverRoute', obj.makeSolverRoute(id), ...
                'plotConfig', obj.makePlotConfig(family), ...
                'resultLayout', obj.makeResultLayout(id));
        end

        function value = makeStatus(obj, id) %#ok<INUSD>
            value = 'available';

            switch id
                case 'SHOCK_OBLIQUE_R'
                    value = 'notIncluded';
            end

        end

        function defaults = makeDefaults(obj, id)
            defaults = obj.blankDefaults();
            defaults.PR1 = '300';
            defaults.PR2 = '1';

            switch id
                case 'TP'
                    defaults.PP1 = '2500';
                    defaults.PP2 = '1';
                case {'HP', 'EV'}
                    defaults.PP2 = '1';
                case 'SP'
                    defaults.PP2 = '100';
                case 'TV'
                    defaults.PP1 = '2500';
                    defaults.PP2 = '1';
                case 'SV'
                    defaults.PP2 = '1';
                    defaults.PP4 = '2';
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR', ...
                        'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_ACOUSTIC'}
                    defaults.PR4 = '2';
                case 'SHOCK_OBLIQUE'
                    defaults.PR4 = '5';
                    defaults.PR5 = '40';
                case 'SHOCK_POLAR_R'
                    defaults.PR1 = '226.65';
                    defaults.PR2 = '0.0117';
                    defaults.PR4 = '20';
                    defaults.PP5 = '35';
                case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R', ...
                        'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R', ...
                        'DET_POLAR'}
                    defaults.PR3 = '2';
                case 'DET_OBLIQUE'
                    defaults.PR3 = '2';
                    defaults.PR4 = '60';
                case 'SHOCKTURBULENCE_VORTICAL_ENTROPIC'
                    defaults.PR4 = '2';
                    defaults.PP1 = '-0.1';
                case 'SHOCKTURBULENCE_COMPRESSIBLE'
                    defaults.PR4 = '2';
                    defaults.PP1 = '0.1';
                    defaults.PP2 = '0';
                    defaults.PP6 = '0.04';
            end

        end

        function defaults = blankDefaults(obj) %#ok<MANU>
            defaults = struct( ...
                'PR1', '', ...
                'PR2', '', ...
                'PR3', '', ...
                'PR4', '', ...
                'PR5', '', ...
                'PP1', '', ...
                'PP2', '', ...
                'PP3', '', ...
                'PP4', '', ...
                'PP5', '', ...
                'PP6', '');
        end

        function labels = makeLabels(obj, id)
            labels = obj.blankLabels();

            switch id
                case 'TV'
                    labels.text_RP2 = 'Specific volume [m3/kg]';
                case 'HP'
                    labels.text_RP3 = 'Constant Enthalpy: hP = hR';
                case 'SP'
                    labels.text_RP3 = 'Constant Entropy: SP = SR';
                case 'EV'
                    labels.text_RP2 = 'Specific volume [m3/kg]';
                    labels.text_RP3 = 'Constant Internal energy: eP = eR';
                    labels.text_RP4 = 'Constant Volume: vP = vR';
                case 'SV'
                    labels.text_RP2 = 'Specific volume [m3/kg]';
                    labels.text_RP3 = 'Constant Entropy: SP = SR';
                    labels.text_RP4 = 'Volume Products/Reactants';
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR', ...
                        'SHOCK_POLAR_R', ...
                        'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_ACOUSTIC'}
                    labels.text_RP3 = 'Shock velocity [m/s]';
                    labels.text_RP4 = 'Mach number [-]';
                case 'SHOCK_OBLIQUE'
                    labels.text_RP = 'Reactants';
                    labels.text_RP3 = 'Shock velocity [m/s]';
                    labels.text_RP4 = 'Mach number [-]';
                    labels.text_RP5 = 'Wave/Deflection angle [deg]';
                case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R'}
                    labels.text_RP3 = 'Overdriven parameter [-]';
                case {'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R'}
                    labels.text_RP3 = 'Underdriven parameter [-]';
                case 'DET_OBLIQUE'
                    labels.text_RP = 'Reactants';
                    labels.text_RP3 = 'Driven parameter [-]';
                    labels.text_RP4 = 'Wave/Deflection angle [deg]';
                case 'DET_POLAR'
                    labels.text_RP3 = 'Overdriven parameter [-]';
                case 'ROCKET'
                    labels.additionalPanelTitle = 'Optional parameters';
                    labels.text_RP3 = 'Area ratios A/A_throat';
                    labels.text_R2 = 'Subsonic';
                    labels.text_P2 = 'Supersonic';
                case 'SHOCKTURBULENCE_VORTICAL_ENTROPIC'
                    labels.text_RP3 = 'Shock velocity [m/s]';
                    labels.text_RP4 = 'Mach number [-]';
                    labels.text_P1 = 'Upstream turbulence';
                    labels.text_RP1_2 = 'chi [-]';
                case 'SHOCKTURBULENCE_COMPRESSIBLE'
                    labels.text_RP3 = 'Shock velocity [m/s]';
                    labels.text_RP4 = 'Mach number [-]';
                    labels.text_P1 = 'Upstream turbulence';
                    labels.text_RP1_2 = 'eta [-]';
                    labels.text_RP2_2 = 'chi [-]';
                    labels.text_RP3_2 = 'etaVorticity [-]';
            end

        end

        function labels = blankLabels(obj) %#ok<MANU>
            labels = struct( ...
                'text_RP2', 'Pressure [bar]', ...
                'additionalPanelTitle', 'Additional constraints', ...
                'text_P1', 'Products', ...
                'text_RP', 'Products', ...
                'text_R2', 'Reactants', ...
                'text_P2', 'Products', ...
                'text_RP3', '', ...
                'text_RP4', '', ...
                'text_RP5', '', ...
                'text_RP1_2', 'Area ratio A_c/A_t', ...
                'text_RP2_2', 'Mass flux [kg/s]', ...
                'text_RP3_2', 'etaVorticity [-]');
        end

        function visibility = makeVisibility(obj, id)
            visibility = obj.blankVisibility();

            switch id
                case 'TP'
                    visibility.PP1 = true;
                    visibility.PP2 = true;
                case {'HP', 'SP'}
                    visibility.PP2 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                case 'TV'
                    visibility.PP1 = true;
                    visibility.PP2 = true;
                case 'EV'
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                    visibility.text_RP4 = true;
                case 'SV'
                    visibility.PP4 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                    visibility.text_RP4 = true;
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR', ...
                        'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_ACOUSTIC'}
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_R2 = true;
                    visibility.text_RP4 = true;
                    visibility.shocks = true;
                case 'SHOCK_OBLIQUE'
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.PR5 = true;
                    visibility.PP5 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                    visibility.text_RP4 = true;
                    visibility.text_RP5 = true;
                    visibility.shocks = true;
                    visibility.oblique = true;
                case 'SHOCK_POLAR_R'
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.PR5 = true;
                    visibility.PP5 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                    visibility.text_RP4 = true;
                    visibility.text_RP5 = true;
                    visibility.shocks = true;
                case {'DET', 'DET_R'}
                    visibility.shocks = true;
                case {'DET_OVERDRIVEN', 'DET_OVERDRIVEN_R', ...
                        'DET_UNDERDRIVEN', 'DET_UNDERDRIVEN_R', ...
                        'DET_POLAR'}
                    visibility.PR3 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_R2 = true;
                    visibility.shocks = true;
                case 'DET_OBLIQUE'
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.PP4 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_RP = true;
                    visibility.text_RP4 = true;
                    visibility.shocks = true;
                    visibility.oblique = true;
                case 'ROCKET'
                    visibility.FLAG_IAC = true;
                    visibility.PR3 = true;
                    visibility.PP3 = true;
                    visibility.text_P1 = false;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_R2 = true;
                    visibility.text_P2 = true;
                    visibility.shocks = true;
                    visibility.rocket = true;
                case 'SHOCKTURBULENCE_VORTICAL_ENTROPIC'
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.PP1 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_P1 = true;
                    visibility.text_R2 = true;
                    visibility.text_RP4 = true;
                    visibility.shocks = true;
                    visibility.shockTurbulence = true;
                    visibility.shockTurbulenceP1 = true;
                    visibility.TurbulenceStatistics = true;
                case 'SHOCKTURBULENCE_COMPRESSIBLE'
                    visibility.PR3 = true;
                    visibility.PR4 = true;
                    visibility.PP1 = true;
                    visibility.PP2 = true;
                    visibility.PP6 = true;
                    visibility.AdditionalconstraintsPanel = true;
                    visibility.text_P1 = true;
                    visibility.text_R2 = true;
                    visibility.text_RP4 = true;
                    visibility.shocks = true;
                    visibility.shockTurbulence = true;
                    visibility.shockTurbulenceP1 = true;
                    visibility.shockTurbulenceP2 = true;
                    visibility.TurbulenceStatistics = true;
            end

            if any(strcmp(id, {'SHOCKTURBULENCE_VORTICAL', 'SHOCKTURBULENCE_ACOUSTIC'}))
                visibility.shockTurbulence = true;
                visibility.TurbulenceStatistics = true;
                visibility.text_P1 = false;
            end

        end

        function visibility = blankVisibility(obj) %#ok<MANU>
            visibility = struct( ...
                'FLAG_IAC', false, ...
                'PP1', false, ...
                'PP2', false, ...
                'PP3', false, ...
                'PP4', false, ...
                'PP5', false, ...
                'PP6', false, ...
                'PR3', false, ...
                'PR4', false, ...
                'PR5', false, ...
                'text_P1', true, ...
                'AdditionalconstraintsPanel', false, ...
                'text_RP', false, ...
                'text_R2', false, ...
                'text_P2', false, ...
                'text_RP4', false, ...
                'text_RP5', false, ...
                'shocks', false, ...
                'oblique', false, ...
                'rocket', false, ...
                'shockTurbulence', false, ...
                'shockTurbulenceP1', false, ...
                'shockTurbulenceP2', false, ...
                'TurbulenceStatistics', false);
        end

        function roles = makeInputRoles(~, id)
            roles = struct( ...
                'reactants', {{'PR1', 'PR2'}}, ...
                'products', {{}}, ...
                'additionalReactants', {{}}, ...
                'additionalProducts', {{}}, ...
                'selector', 'none');

            switch id
                case {'TP', 'TV'}
                    roles.products = {'PP1', 'PP2'};
                case {'HP', 'SP', 'EV'}
                    roles.products = {'PP2'};
                case 'SV'
                    roles.products = {'PP4'};
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_POLAR', ...
                        'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_ACOUSTIC'}
                    roles.additionalReactants = {'PR3', 'PR4'};
                    roles.selector = 'mach';
                case {'SHOCK_OBLIQUE', 'SHOCK_POLAR_R'}
                    roles.additionalReactants = {'PR3', 'PR4', 'PR5'};
                    roles.additionalProducts = {'PP5'};
                    roles.selector = 'betaOrTheta';
                case {'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', ...
                        'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R', ...
                        'DET_POLAR'}
                    roles.additionalReactants = {'PR3'};
                    roles.selector = 'driveFactor';
                case 'DET_OBLIQUE'
                    roles.additionalReactants = {'PR3', 'PR4'};
                    roles.additionalProducts = {'PP4'};
                    roles.selector = 'driveFactorBetaOrTheta';
                case 'ROCKET'
                    roles.additionalReactants = {'PR3'};
                    roles.additionalProducts = {'PP1', 'PP3'};
                    roles.selector = 'rocketModel';
                case 'SHOCKTURBULENCE_VORTICAL_ENTROPIC'
                    roles.additionalReactants = {'PR3', 'PR4'};
                    roles.additionalProducts = {'PP1'};
                    roles.selector = 'machChi';
                case 'SHOCKTURBULENCE_COMPRESSIBLE'
                    roles.additionalReactants = {'PR3', 'PR4'};
                    roles.additionalProducts = {'PP1', 'PP2', 'PP6'};
                    roles.selector = 'machEtaChiEtaVorticity';
            end

        end

        function route = makeSolverRoute(obj, id) %#ok<INUSD>
            route = struct( ...
                'solver', 'equilibriumSolver', ...
                'problemType', id, ...
                'solverProblemType', id, ...
                'variantPolicy', 'none', ...
                'outputMixtures', 2, ...
                'maxOutputMixtures', 2, ...
                'requiresSubsolverSilence', false);

            switch id
                case {'TP', 'HP', 'SP', 'TV', 'EV', 'SV'}
                    route.solver = 'equilibriumSolver';
                case {'SHOCK_I', 'SHOCK_R', 'SHOCK_OBLIQUE', ...
                        'SHOCK_OBLIQUE_R', 'SHOCK_POLAR', 'SHOCK_POLAR_R'}
                    route.solver = 'shockSolver';
                    route.requiresSubsolverSilence = true;
                case {'DET', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', ...
                        'DET_R', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R', ...
                        'DET_OBLIQUE', 'DET_POLAR'}
                    route.solver = 'detonationSolver';
                    route.requiresSubsolverSilence = true;
                case 'ROCKET'
                    route.solver = 'rocketSolver';
                    route.solverProblemType = 'ROCKET_IAC';
                    route.variantPolicy = 'rocketModelAndAreaRatio';
                    route.outputMixtures = 3;
                    route.maxOutputMixtures = 5;
                    route.requiresSubsolverSilence = true;
                case {'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_VORTICAL_ENTROPIC', ...
                        'SHOCKTURBULENCE_ACOUSTIC', ...
                        'SHOCKTURBULENCE_COMPRESSIBLE'}
                    route.solver = 'shockTurbulenceSolver';
                    route.solverProblemType = strrep(id, 'SHOCKTURBULENCE_', '');
                    route.outputMixtures = 2;
                    route.requiresSubsolverSilence = true;
            end

            switch id
                case 'SHOCK_OBLIQUE'
                    route.variantPolicy = 'betaOrTheta';
                    route.outputMixtures = 2;
                    route.maxOutputMixtures = 3;
                case 'SHOCK_OBLIQUE_R'
                    route.variantPolicy = 'betaOrTheta';
                    route.outputMixtures = 4;
                    route.maxOutputMixtures = 4;
                case 'SHOCK_POLAR_R'
                    route.variantPolicy = 'betaOrTheta';
                    route.outputMixtures = 6;
                    route.maxOutputMixtures = 6;
                case 'DET_OBLIQUE'
                    route.variantPolicy = 'betaOrTheta';
                    route.outputMixtures = 3;
                    route.maxOutputMixtures = 3;
                case {'SHOCK_R', 'DET_R', 'DET_OVERDRIVEN_R', ...
                        'DET_UNDERDRIVEN_R'}
                    route.outputMixtures = 3;
                    route.maxOutputMixtures = 3;
            end

        end

        function config = makePlotConfig(obj, family) %#ok<INUSD>
            config = struct( ...
                'properties', {{'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}}, ...
                'basis', {{[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}});

            switch family
                case 'detonation'
                    config.properties = {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'uShock'};
                    config.basis = {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], []};
                case 'rocket'
                    config.properties = {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'};
                    config.basis = {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []};
                case 'shockturbulence'
                    config.properties = {'K', 'R11', 'RTT', 'Ka', 'Kr', 'enstrophy', 'kolmogorovLengthRatio'};
                    config.basis = {[], [], [], [], [], [], []};
            end

        end

        function layout = makeResultLayout(obj, id)
            layout = struct( ...
                'mode', 'reactantsProducts', ...
                'states', {{'mix1', 'mix2'}}, ...
                'selectionState', 'mix2', ...
                'variantStates', struct(), ...
                'variantOutputStates', struct());
            layout.outputStates = obj.makeOutputStates( ...
                {'reactants', 'products'}, ...
                {'Reactants', 'Products'}, ...
                {'mix1', 'mix2'}, ...
                {'initial', 'equilibrium'});

            switch id
                case 'SHOCK_I'
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preShock', 'postShock'}, ...
                        {'Pre-shock', 'Post-shock'}, ...
                        {'mix1', 'mix2'}, ...
                        {'initial', 'branch'});
                case {'DET_R', ...
                        'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R'}
                    layout.states = {'mix1', 'mix2', 'mix3'};
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preDetonation', 'postDetonation', 'reflectedPostDetonation'}, ...
                        {'Pre-detonation', 'Incident post-detonation', 'Reflected post-detonation'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'branch', 'branch'});
                case 'SHOCK_R'
                    layout.states = {'mix1', 'mix2', 'mix3'};
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preShock', 'postShock', 'reflectedPostShock'}, ...
                        {'Pre-shock', 'Incident post-shock', 'Reflected post-shock'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'branch', 'branch'});
                case 'SHOCK_OBLIQUE'
                    layout.variantStates.beta = {'mix1', 'mix2'};
                    layout.variantStates.theta = {'mix1', 'mix2', 'mix3'};
                    layout.variantOutputStates.beta = obj.makeOutputStates( ...
                        {'preShock', 'postShock'}, ...
                        {'Pre-shock', 'Post-shock'}, ...
                        {'mix1', 'mix2'}, ...
                        {'initial', 'branch'});
                    layout.variantOutputStates.theta = obj.makeOutputStates( ...
                        {'preShock', 'weakPostShock', 'strongPostShock'}, ...
                        {'Pre-shock', 'Weak post-shock', 'Strong post-shock'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'branch', 'branch'});
                    layout.outputStates = layout.variantOutputStates.beta;
                case 'SHOCK_OBLIQUE_R'
                    layout.states = {'mix1', 'mix2', 'mix3', 'mix4'};
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preShock', 'postShock', 'weakReflectedPostShock', 'strongReflectedPostShock'}, ...
                        {'Pre-shock', 'Incident post-shock', 'Weak reflected post-shock', 'Strong reflected post-shock'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4'}, ...
                        {'initial', 'branch', 'branch', 'branch'});
                case 'SHOCK_POLAR'
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preShock', 'postShockPolar'}, ...
                        {'Pre-shock', 'Post-shock polar'}, ...
                        {'mix1', 'mix2'}, ...
                        {'initial', 'polar'});
                case 'SHOCK_POLAR_R'
                    layout.states = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5', 'mix6'};
                    layout.variantStates.beta = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5'};
                    layout.variantStates.theta = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5', 'mix6'};
                    layout.variantOutputStates.beta = obj.makeOutputStates( ...
                        {'preShock', 'incidentPolar', 'postShock', 'reflectedPolar', 'reflectedPostShock'}, ...
                        {'Pre-shock', 'Incident polar', 'Incident post-shock', 'Reflected polar', 'Reflected post-shock'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4', 'mix5'}, ...
                        {'initial', 'polar', 'branch', 'polar', 'branch'});
                    layout.variantOutputStates.theta = obj.makeOutputStates( ...
                        {'preShock', 'incidentPolar', 'weakPostShock', 'reflectedPolar', 'weakReflectedPostShock', 'strongReflectedPostShock'}, ...
                        {'Pre-shock', 'Incident polar', 'Weak incident post-shock', 'Reflected polar', 'Weak reflected post-shock', 'Strong reflected post-shock'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4', 'mix5', 'mix6'}, ...
                        {'initial', 'polar', 'branch', 'polar', 'branch', 'branch'});
                    layout.outputStates = layout.variantOutputStates.beta;
                case {'DET', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN'}
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preDetonation', 'postDetonation'}, ...
                        {'Pre-detonation', 'Post-detonation'}, ...
                        {'mix1', 'mix2'}, ...
                        {'initial', 'branch'});
                case 'DET_OBLIQUE'
                    layout.states = {'mix1', 'mix2', 'mix3'};
                    layout.variantOutputStates.beta = obj.makeOutputStates( ...
                        {'preDetonation', 'underdrivenPostDetonation', 'overdrivenPostDetonation'}, ...
                        {'Pre-detonation', 'Underdriven post-detonation', 'Overdriven post-detonation'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'branch', 'branch'});
                    layout.variantOutputStates.theta = obj.makeOutputStates( ...
                        {'preDetonation', 'weakPostDetonation', 'strongPostDetonation'}, ...
                        {'Pre-detonation', 'Weak post-detonation', 'Strong post-detonation'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'branch', 'branch'});
                    layout.outputStates = layout.variantOutputStates.beta;
                case 'DET_POLAR'
                    layout.outputStates = obj.makeOutputStates( ...
                        {'preDetonation', 'postDetonationPolar'}, ...
                        {'Pre-detonation', 'Post-detonation polar'}, ...
                        {'mix1', 'mix2'}, ...
                        {'initial', 'polar'});
                case 'ROCKET'
                    layout.mode = 'rocket';
                    layout.states = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5'};
                    layout.variantStates.IAC = {'mix1', 'mix2', 'mix3'};
                    layout.variantStates.IAC_ARATIO = {'mix1', 'mix2', 'mix3', 'mix4'};
                    layout.variantStates.FAC = {'mix1', 'mix2', 'mix3', 'mix4'};
                    layout.variantStates.FAC_ARATIO = {'mix1', 'mix2', 'mix3', 'mix4', 'mix5'};
                    layout.variantOutputStates.IAC = obj.makeOutputStates( ...
                        {'reactants', 'chamber', 'throat'}, ...
                        {'Reactants', 'Chamber', 'Throat'}, ...
                        {'mix1', 'mix2', 'mix3'}, ...
                        {'initial', 'station', 'station'});
                    layout.variantOutputStates.IAC_ARATIO = obj.makeOutputStates( ...
                        {'reactants', 'chamber', 'throat', 'exit'}, ...
                        {'Reactants', 'Chamber', 'Throat', 'Exit'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4'}, ...
                        {'initial', 'station', 'station', 'station'});
                    layout.variantOutputStates.FAC = obj.makeOutputStates( ...
                        {'reactants', 'injector', 'chamber', 'throat'}, ...
                        {'Reactants', 'Injector', 'Chamber', 'Throat'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4'}, ...
                        {'initial', 'station', 'station', 'station'});
                    layout.variantOutputStates.FAC_ARATIO = obj.makeOutputStates( ...
                        {'reactants', 'injector', 'chamber', 'throat', 'exit'}, ...
                        {'Reactants', 'Injector', 'Chamber', 'Throat', 'Exit'}, ...
                        {'mix1', 'mix2', 'mix3', 'mix4', 'mix5'}, ...
                        {'initial', 'station', 'station', 'station', 'station'});
                    layout.outputStates = layout.variantOutputStates.IAC;
                case {'SHOCKTURBULENCE_VORTICAL', ...
                        'SHOCKTURBULENCE_VORTICAL_ENTROPIC', ...
                        'SHOCKTURBULENCE_ACOUSTIC', ...
                        'SHOCKTURBULENCE_COMPRESSIBLE'}
                    layout.mode = 'shockTurbulence';
                    layout.states = {'resultsLIA', 'mix1', 'mix2'};
                    layout.outputStates = obj.makeOutputStates( ...
                        {'statistics', 'preShock', 'postShock'}, ...
                        {'Turbulence statistics', 'Pre-shock', 'Post-shock'}, ...
                        {'resultsLIA', 'mix1', 'mix2'}, ...
                        {'statistics', 'initial', 'branch'}, ...
                        {'data', 'mixture', 'mixture'});
            end

        end

        function outputStates = makeOutputStates(~, ids, labels, fields, roles, types)
            if nargin < 5 || isempty(roles)
                roles = ids;
            end

            if nargin < 6 || isempty(types)
                types = repmat({'mixture'}, size(ids));
            end

            outputStates = repmat(struct( ...
                'id', '', ...
                'label', '', ...
                'field', '', ...
                'role', '', ...
                'type', ''), 1, numel(ids));

            for i = 1:numel(ids)
                outputStates(i).id = ids{i};
                outputStates(i).label = labels{i};
                outputStates(i).field = fields{i};
                outputStates(i).role = roles{i};
                outputStates(i).type = types{i};
            end
        end

        function key = normalize(obj, problemType) %#ok<INUSD>
            key = matlab.lang.makeValidName(char(problemType));
        end
    end
end

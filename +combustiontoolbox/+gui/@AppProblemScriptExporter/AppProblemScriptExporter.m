classdef AppProblemScriptExporter < handle
    % Converts canonical GUI problem setups into runnable MATLAB scripts.
    %
    % Attributes:
    %     problemCatalog (AppProblemCatalog): Problem metadata and outputs
    %
    % Examples:
    %     * exporter = combustiontoolbox.gui.AppProblemScriptExporter();
    %     * script = exporter.export(problem.setup);

    properties (Access = private)
        problemCatalog
    end

    methods
        function obj = AppProblemScriptExporter(varargin)
            % AppProblemScriptExporter constructor
            %
            % Optional Args:
            %     problemCatalog (AppProblemCatalog): Problem metadata catalog
            %
            % Returns:
            %     obj (AppProblemScriptExporter): Initialized script exporter
            obj.problemCatalog = combustiontoolbox.gui.AppProblemCatalog();

            if nargin > 0 && isa(varargin{1}, 'combustiontoolbox.gui.AppProblemCatalog')
                obj.problemCatalog = varargin{1};
            end
        end

        function script = export(obj, setup)
            % Export a canonical problem setup to a MATLAB script
            %
            % Args:
            %     setup (struct): Canonical problem setup
            %
            % Returns:
            %     script (char): Runnable MATLAB script text
            obj.validateSetup(setup);

            lines = [ ...
                obj.headerLines(setup), ...
                obj.importLines(setup), ...
                obj.databaseLines(), ...
                obj.chemicalSystemLines(setup), ...
                obj.mixtureLines(), ...
                obj.compositionLines(setup), ...
                obj.propertiesLines(setup), ...
                obj.solverLines(setup), ...
                obj.postprocessLines(setup)];

            script = strjoin(lines, newline);
            script = [script, newline];
        end
    end

    methods (Access = private)
        function validateSetup(~, setup)
            if ~isstruct(setup)
                error('AppProblemScriptExporter:InvalidSetup', ...
                    'Script export requires a canonical problem setup struct.');
            end

            requiredFields = {'problemType', 'requestedProblemType', 'family', ...
                'reactants', 'productSpecies', 'reactantStateProperties', ...
                'productStateProperties', 'additionalInputsReactants', ...
                'additionalInputsProducts', 'options'};

            for i = 1:numel(requiredFields)
                if ~isfield(setup, requiredFields{i})
                    error('AppProblemScriptExporter:InvalidSetup', ...
                        'Script export setup is missing the "%s" field.', requiredFields{i});
                end
            end
        end

        function lines = headerLines(~, setup)
            lines = { ...
                '% -------------------------------------------------------------------------', ...
                ['% EXAMPLE: ', char(setup.problemType)], ...
                '%', ...
                '% Generated from Combustion Toolbox App', ...
                '% -------------------------------------------------------------------------', ...
                ''};
        end

        function lines = importLines(obj, setup)
            lines = { ...
                '% Import packages', ...
                'import combustiontoolbox.databases.NasaDatabase', ...
                'import combustiontoolbox.core.*', ...
                ['import ', obj.solverPackage(setup), '.*'], ...
                ''};
        end

        function lines = databaseLines(~)
            lines = { ...
                '% Get database', ...
                'DB = NasaDatabase();', ...
                ''};
        end

        function lines = chemicalSystemLines(obj, setup)
            lines = { ...
                '% Set chemical system', ...
                sprintf('system = ChemicalSystem(%s);', ...
                    obj.chemicalSystemArguments(setup)), ...
                ''};
        end

        function text = chemicalSystemArguments(obj, setup)
            tokens = {'DB'};
            optionTokens = obj.chemicalSystemOptionArguments(setup);

            if ~obj.usesAutomaticProductSpecies(setup)
                tokens{end + 1} = obj.scriptValue(setup.productSpecies);
            elseif ~isempty(optionTokens)
                tokens{end + 1} = '[]';
            end

            tokens = [tokens, optionTokens];
            text = strjoin(tokens, ', ');
        end

        function tokens = chemicalSystemOptionArguments(obj, setup)
            tokens = {};

            if obj.setupFlag(setup, 'ionizedSpecies', false)
                tokens = [tokens, {obj.scriptValue('FLAG_ION'), ...
                    obj.scriptLogical(true)}]; %#ok<AGROW>
            end
        end

        function value = usesAutomaticProductSpecies(~, setup)
            value = false;

            if ~isfield(setup, 'products') || ~isfield(setup.products, 'type')
                return
            end

            value = strcmpi(setup.products.type, 'auto');
        end

        function value = setupFlag(~, setup, name, defaultValue)
            value = defaultValue;

            if isfield(setup, 'flags') && isfield(setup.flags, name)
                value = setup.flags.(name);
            end
        end

        function lines = mixtureLines(~)
            lines = { ...
                '% Initialize mixture', ...
                'mix = Mixture(system);', ...
                ''};
        end

        function lines = compositionLines(obj, setup)
            lines = {'% Set chemical state'};

            if isempty(setup.reactants.species)
                lines{end + 1} = '';
                return
            end

            species = setup.reactants.species;
            amounts = setup.reactants.amount;
            hasMixedTemperatures = obj.hasMixedReactantSpeciesTemperatures(setup);

            if isempty(amounts)
                amounts = ones(1, numel(species));
            end

            if isempty(setup.reactants.role)
                lines{end + 1} = sprintf('set(mix, %s, %s);', ...
                    obj.scriptValue(species), obj.scriptValue(amounts));

                if hasMixedTemperatures
                    lines{end + 1} = sprintf('setTemperatureSpecies(mix, %s);', ...
                        obj.scriptValue(setup.reactants.temperature));
                end
            else
                if hasMixedTemperatures
                    temperatures = setup.reactants.temperature;
                else
                    temperatures = [];
                end

                roleLines = obj.compositionRoleLines(setup, species, amounts, temperatures);
                lines = [lines, roleLines];
            end

            lines{end + 1} = '';
        end

        function lines = compositionRoleLines(obj, setup, species, amounts, temperatures)
            roles = setup.reactants.role;
            roleNames = {'fuel', 'oxidizer', 'inert'};
            assigned = false(size(species));
            hasMixedTemperatures = ~isempty(temperatures);
            lines = cell(1, numel(roleNames) + 2);
            temperatureList = zeros(1, numel(species));
            count = 0;
            temperatureCount = 0;

            for i = 1:numel(roleNames)
                index = strcmpi(roles, roleNames{i});

                if ~any(index)
                    continue
                end

                assigned(index) = true;
                count = count + 1;
                lines{count} = sprintf('set(mix, %s, %s, %s);', ...
                    obj.scriptValue(species(index)), ...
                    obj.scriptValue(roleNames{i}), ...
                    obj.scriptValue(amounts(index)));

                if hasMixedTemperatures
                    nextTemperatureCount = temperatureCount + nnz(index);
                    temperatureList(temperatureCount + 1:nextTemperatureCount) = temperatures(index);
                    temperatureCount = nextTemperatureCount;
                end
            end

            index = ~assigned;

            if any(index)
                count = count + 1;
                lines{count} = sprintf('set(mix, %s, %s);', ...
                    obj.scriptValue(species(index)), obj.scriptValue(amounts(index)));

                if hasMixedTemperatures
                    nextTemperatureCount = temperatureCount + nnz(index);
                    temperatureList(temperatureCount + 1:nextTemperatureCount) = temperatures(index);
                    temperatureCount = nextTemperatureCount;
                end
            end

            if hasMixedTemperatures
                temperatureList = temperatureList(1:temperatureCount);
                count = count + 1;
                lines{count} = sprintf('setTemperatureSpecies(mix, %s);', ...
                    obj.scriptValue(temperatureList));
            end

            lines = lines(1:count);
        end

        function value = hasMixedReactantSpeciesTemperatures(~, setup)
            value = isfield(setup.reactants, 'temperature') ...
                && numel(setup.reactants.temperature) > 1 ...
                && any(abs(setup.reactants.temperature - setup.reactants.temperature(1)) > 1e-12);
        end

        function lines = propertiesLines(obj, setup)
            if obj.hasMixedReactantSpeciesTemperatures(setup)
                ignoredProperties = {'temperature', 't'};
            else
                ignoredProperties = {};
            end

            reactantArguments = obj.setupArguments( ...
                setup.reactantStateProperties, setup.additionalInputsReactants, ...
                ignoredProperties, setup);
            lines = { ...
                '% Set properties', ...
                sprintf('mixArrayReactants = setProperties(mix, %s);', reactantArguments)};

            if obj.isEquilibrium(setup)
                productArguments = obj.setupArguments( ...
                    setup.productStateProperties, setup.additionalInputsProducts, ...
                    {}, setup);
                lines{end + 1} = sprintf('mixArrayProducts = setProperties(mix, %s);', productArguments);
            end

            lines{end + 1} = '';
        end

        function lines = solverLines(obj, setup)
            solverClass = obj.solverClass(setup);
            solverProblemType = obj.solverProblemType(setup);
            printResults = obj.scriptLogical(setup.options.printResults);
            lines = { ...
                '% Set problem', ...
                sprintf('solver = %s(''problemType'', %s, ''FLAG_RESULTS'', %s);', ...
                    solverClass, obj.scriptValue(solverProblemType), printResults)};

            if obj.isRocket(setup) && isfield(setup.flags, 'rocketSubsonic') ...
                    && ~isempty(setup.flags.rocketSubsonic)
                lines{end + 1} = sprintf('solver.FLAG_SUBSONIC = %s;', ...
                    obj.scriptLogical(setup.flags.rocketSubsonic));
            end

            lines{end + 1} = '';
            lines{end + 1} = '% Solve problem';

            if obj.isEquilibrium(setup)
                lines{end + 1} = 'solver.solveArray(mixArrayProducts);';
            else
                outputVariables = obj.outputVariables(setup);
                solverMethod = 'solver.solveArray';

                if obj.isShockTurbulence(setup)
                    solverMethod = 'solver.solve';
                end

                lines{end + 1} = obj.solverCallLine( ...
                    solverMethod, 'mixArrayReactants', outputVariables);
            end

            lines{end + 1} = '';
        end

        function lines = postprocessLines(obj, setup)
            lines = {};

            if ~isfield(setup.options, 'reportType') || strcmpi(setup.options.reportType, 'None')
                return
            end

            if obj.isEquilibrium(setup)
                reportArguments = {'mixArrayProducts'};
            else
                reportArguments = obj.outputVariables(setup);
            end

            lines = { ...
                '% Generate report', ...
                sprintf('report(solver, %s);', strjoin(reportArguments, ', '))};
        end

        function text = setupArguments(obj, stateProperties, additionalInputs, varargin)
            ignoredProperties = {};
            setup = struct();

            if nargin > 3
                ignoredProperties = varargin{1};
            end

            if nargin > 4
                setup = varargin{2};
            end

            tokens = cell(1, 2 * numel(stateProperties) + numel(additionalInputs));
            count = 0;

            for i = 1:numel(stateProperties)
                if any(strcmpi(stateProperties(i).name, ignoredProperties))
                    continue
                end

                if isempty(stateProperties(i).value)
                    continue
                end

                tokens{count + 1} = obj.scriptValue(stateProperties(i).name);
                tokens{count + 2} = obj.scriptInputValue(setup, stateProperties(i).value);
                count = count + 2;
            end

            for i = 1:2:numel(additionalInputs)
                if i + 1 > numel(additionalInputs) || isempty(additionalInputs{i + 1})
                    continue
                end

                tokens{count + 1} = obj.scriptValue(additionalInputs{i});
                tokens{count + 2} = obj.scriptInputValue(setup, additionalInputs{i + 1});
                count = count + 2;
            end

            text = strjoin(tokens(1:count), ', ');
        end

        function text = scriptInputValue(obj, setup, value)
            if isfield(setup, 'expressions')
                for i = 1:numel(setup.expressions)
                    expression = setup.expressions(i);

                    if isequal(size(expression.value), size(value)) ...
                            && all(abs(expression.value(:) - value(:)) < 1e-12)
                        text = expression.text;
                        return
                    end
                end
            end

            text = obj.scriptValue(value);
        end

        function line = solverCallLine(~, solverMethod, inputVariable, outputVariables)
            if isscalar(outputVariables)
                line = sprintf('%s = %s(%s);', ...
                    outputVariables{1}, solverMethod, inputVariable);
            else
                line = sprintf('[%s] = %s(%s);', ...
                    strjoin(outputVariables, ', '), solverMethod, inputVariable);
            end
        end

        function outputVariables = outputVariables(obj, setup)
            states = obj.outputStates(setup);
            outputVariables = cell(1, numel(states));

            for i = 1:numel(states)
                outputVariables{i} = matlab.lang.makeValidName(states(i).id);
            end
        end

        function states = outputStates(obj, setup)
            layout = obj.problemCatalog.resultLayout(setup.requestedProblemType);
            states = layout.outputStates;

            if ~isfield(layout, 'variantOutputStates')
                return
            end

            variant = char(setup.problemType);
            variantNames = fieldnames(layout.variantOutputStates);

            for i = 1:numel(variantNames)
                variantName = variantNames{i};

                if strcmpi(variant, variantName) || endsWith(upper(variant), ['_', upper(variantName)])
                    states = layout.variantOutputStates.(variantName);
                    return
                end
            end
        end

        function value = solverPackage(obj, setup)
            switch obj.family(setup)
                case 'equilibrium'
                    value = 'combustiontoolbox.equilibrium';
                case {'shock', 'detonation'}
                    value = 'combustiontoolbox.shockdetonation';
                case 'rocket'
                    value = 'combustiontoolbox.rocket';
                case 'shockturbulence'
                    value = 'combustiontoolbox.shockturbulence';
                otherwise
                    error('AppProblemScriptExporter:UnsupportedFamily', ...
                        'Script export does not support the "%s" problem family.', setup.family);
            end
        end

        function value = solverClass(obj, setup)
            switch obj.family(setup)
                case 'equilibrium'
                    value = 'EquilibriumSolver';
                case 'shock'
                    value = 'ShockSolver';
                case 'detonation'
                    value = 'DetonationSolver';
                case 'rocket'
                    value = 'RocketSolver';
                case 'shockturbulence'
                    value = 'ShockTurbulenceSolver';
                otherwise
                    error('AppProblemScriptExporter:UnsupportedFamily', ...
                        'Script export does not support the "%s" problem family.', setup.family);
            end
        end

        function value = solverProblemType(obj, setup)
            value = char(setup.problemType);
            value = strrep(value, '_BETA', '');
            value = strrep(value, '_THETA', '');
            value = strrep(value, '_ARATIO', '');

            if obj.isShockTurbulence(setup)
                value = lower(strrep(value, 'SHOCKTURBULENCE_', ''));
            end
        end

        function text = scriptValue(obj, value)
            if isempty(value)
                text = '[]';
            elseif ischar(value)
                text = obj.scriptChar(value);
            elseif isstring(value)
                if isscalar(value)
                    text = obj.scriptChar(char(value));
                else
                    text = obj.scriptCellString(cellstr(value));
                end
            elseif iscell(value)
                text = obj.scriptCellString(value);
            elseif islogical(value)
                if isscalar(value)
                    text = obj.scriptLogical(value);
                else
                    text = mat2str(value);
                end
            elseif isnumeric(value)
                if isscalar(value)
                    text = sprintf('%.15g', value);
                else
                    text = mat2str(value, 15);
                end
            else
                error('AppProblemScriptExporter:UnsupportedValue', ...
                    'Script export does not support values of class "%s".', class(value));
            end
        end

        function text = scriptCellString(obj, values)
            if isempty(values)
                text = '{}';
                return
            end

            items = cell(1, numel(values));

            for i = 1:numel(values)
                items{i} = obj.scriptValue(values{i});
            end

            text = ['{', strjoin(items, ', '), '}'];
        end

        function text = scriptChar(~, value)
            text = ['''', strrep(value, '''', ''''''), ''''];
        end

        function text = scriptLogical(~, value)
            if value
                text = 'true';
            else
                text = 'false';
            end
        end

        function family = family(~, setup)
            family = lower(char(setup.family));
        end

        function value = isEquilibrium(obj, setup)
            value = strcmp(obj.family(setup), 'equilibrium');
        end

        function value = isRocket(obj, setup)
            value = strcmp(obj.family(setup), 'rocket');
        end

        function value = isShockTurbulence(obj, setup)
            value = strcmp(obj.family(setup), 'shockturbulence');
        end
    end
end

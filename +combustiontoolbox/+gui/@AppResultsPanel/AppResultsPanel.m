classdef AppResultsPanel < handle
    % Projects solved results onto GUI trees, fields, tables, plots, and exports
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %
    % Examples:
    %     * panel = combustiontoolbox.gui.AppResultsPanel(app);
    %     * panel.applySolution(solution);

    properties (Access = private)
        app
        currentResults = struct([])
    end

    methods
        function obj = AppResultsPanel(app)
            % AppResultsPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppResultsPanel): Initialized results panel updater
            if nargin
                obj.app = app;
            end
        end

        function applySolution(obj, solution)
            % Apply a solved study to the results controls
            %
            % Args:
            %     solution (struct): Solver output and GUI result metadata
            if isempty(obj.app) || isempty(solution.results)
                return
            end

            results = solution.results;
            obj.currentResults = results;
            obj.addResultNodes(results);
            obj.updateResult(results(1), false);
            obj.updateCustomFigures(results);
        end

        function results = currentResultsData(obj)
            % Return current GUI result data
            %
            % Returns:
            %     results (struct): Current GUI result data
            results = obj.currentResults;
        end

        function onTreeSelectionChanged(obj)
            % Apply the selected results tree node to the GUI
            if isempty(obj.app) || isempty(obj.app.Tree.SelectedNodes)
                return
            end

            selectedNode = obj.app.Tree.SelectedNodes(1);

            if ~isempty(selectedNode.Children) || isempty(selectedNode.NodeData)
                return
            end

            result = selectedNode.NodeData;

            if ~isstruct(result) || ~isfield(result, 'outputStates')
                return
            end

            obj.updateResult(result, true);
            obj.writeSelectedCaseReport(result, selectedNode.Text);
        end

        function plotCustomFigure(obj)
            % Plot the selected variables in the custom figures tab
            import combustiontoolbox.utils.display.*
            import combustiontoolbox.utils.extensions.brewermap

            if isempty(obj.app) || isempty(obj.currentResults)
                return
            end

            cla(obj.app.UIAxes);
            mixtureFields = obj.checkedNodeText('Tree_mixtures', 'Mixtures');
            xFields = obj.checkedNodeText('Tree_variable_x');
            yFields = obj.checkedNodeText('Tree_variable_y');

            if isempty(mixtureFields) || isempty(xFields) || isempty(yFields)
                return
            end

            xField = xFields{1};
            config = obj.currentPlotConfig(xField, yFields);
            lineStyles = {'-', '--', ':', '-.'};
            colorPalette = config.colorline;

            if numel(yFields) > 1
                colorPalette = brewermap(numel(yFields), config.colorpalette);
            end

            if obj.defaultPlotSettingsEnabled()
                setFigure(obj.app.UIAxes, config);
            end

            legendNames = cell(1, numel(mixtureFields) * numel(yFields));
            legendIndex = numel(legendNames);
            results = obj.currentResults;

            for i = numel(mixtureFields):-1:1
                mixtures = obj.mixtureSeries(results, mixtureFields{i});

                if isempty(mixtures)
                    continue
                end

                xValues = cell2vector(mixtures, xField);
                basis = cell2vector(mixtures, 'mi');
                xValues = obj.convertPlotUnits(xField, xValues, basis);

                if obj.isCompositionProperty(yFields{1})
                    mixtureArray = [mixtures{:}];
                    plotComposition(mixtureArray(1), mixtureArray, xField, yFields{1}, ...
                        'config', config, 'y_var', mixtureArray, 'axes', obj.app.UIAxes);
                    return
                end

                for j = numel(yFields):-1:1
                    yValues = cell2vector(mixtures, yFields{j});
                    yValues = obj.convertPlotUnits(yFields{j}, yValues, basis);
                    color = colorPalette(min(j, size(colorPalette, 1)), :);
                    lineStyle = lineStyles{mod(i - 1, numel(lineStyles)) + 1};
                    plotFigure(xField, xValues, yFields{j}, yValues, ...
                        'config', config, 'axes', obj.app.UIAxes, ...
                        'linestyle', lineStyle, 'color', color);
                    legendNames{legendIndex} = [mixtureFields{i}, ' - ', yFields{j}];
                    legendIndex = legendIndex - 1;
                end
            end

            obj.applyLegend(legendNames, numel(mixtureFields), numel(yFields), config);
        end

        function exportResults(obj, format)
            % Export current mixture output states
            if isempty(obj.app) || isempty(obj.currentResults)
                return
            end

            mixArrays = obj.mixtureOutputArrays(obj.currentResults);

            if isempty(mixArrays) || ~isprop(obj.app, 'export') || isempty(obj.app.export)
                return
            end

            previousFormat = obj.app.export.format;

            try
                obj.app.export.format = format;
                obj.app.export.export(mixArrays{:});
            catch exception
                obj.app.export.format = previousFormat;
                rethrow(exception);
            end

            obj.app.export.format = previousFormat;
        end

        function clear(obj)
            % Clear result controls owned by the results panel
            if isempty(obj.app)
                return
            end

            obj.clearCurrentResults();
            obj.clearAxes();
            obj.clearTreeChildren({'Node_Results', 'Mixtures', 'Variable_x', 'Variable_y'});
            obj.clearCheckedNodes({'Tree_mixtures', 'Tree_variable_x', 'Tree_variable_y'});
            obj.clearProductsTable();
            obj.clearNumericFields();
        end

        function value = problemErrorValue(obj)
            % Return the current problem error value
            %
            % Returns:
            %     value (double): Problem relative error value
            value = 0;

            if isempty(obj.app) || ~isprop(obj.app, 'text_error_problem') ...
                    || isempty(obj.app.text_error_problem)
                return
            end

            value = obj.app.text_error_problem.Value;
        end
    end

    methods (Access = private)
        function clearCurrentResults(obj)
            obj.currentResults = struct([]);
        end

        function clearAxes(obj)
            if isprop(obj.app, 'UIAxes') && ~isempty(obj.app.UIAxes) && isvalid(obj.app.UIAxes)
                cla(obj.app.UIAxes);
            end
        end

        function clearTreeChildren(obj, componentNames)
            for i = 1:numel(componentNames)
                componentName = componentNames{i};

                if isprop(obj.app, componentName) && ~isempty(obj.app.(componentName))
                    delete(obj.app.(componentName).Children);
                end
            end
        end

        function clearCheckedNodes(obj, componentNames)
            for i = 1:numel(componentNames)
                componentName = componentNames{i};

                if isprop(obj.app, componentName) && isprop(obj.app.(componentName), 'CheckedNodes')
                    obj.app.(componentName).CheckedNodes = [];
                end
            end
        end

        function clearProductsTable(obj)
            if isprop(obj.app, 'UITable_P')
                obj.app.UITable_P.Data = {};
            end
        end

        function clearNumericFields(obj)
            prefixes = {'text_T', 'text_p', 'text_r', 'text_h', 'text_e', ...
                'text_cp', 'text_s', 'text_gamma', 'text_W', 'text_sound', ...
                'text_u', 'text_M', 'text_Aratio', 'text_Cstar', 'text_Ivac', ...
                'text_Isp', 'text_beta_min', 'text_beta', 'text_theta'};

            for i = 1:numel(prefixes)
                for suffix = 1:7
                    obj.setNumericValue([prefixes{i}, '_', num2str(suffix)], 0);
                end
            end

            obj.setNumericValue('text_error_problem', 0);
        end

        function addResultNodes(obj, results)
            problemNode = obj.findOrAddNode(obj.app.Node_Results, ...
                ['Problem Type: ', results(1).ProblemType]);
            reactantsNode = obj.findOrAddNode(problemNode, ...
                ['Reactants: ', results(1).Reactants]);
            productsNode = obj.findOrAddNode(reactantsNode, ...
                ['List Products: ', results(1).Products]);

            for i = 1:numel(results)
                label = sprintf('Case %d', productsNode.UserData);
                uitreenode(productsNode, 'Text', label, 'NodeData', results(i));
                productsNode.UserData = productsNode.UserData + 1;
            end

            expand(obj.app.Node_Results, 'all');
        end

        function node = findOrAddNode(~, parentNode, text)
            for i = 1:numel(parentNode.Children)
                if strcmpi(parentNode.Children(i).Text, text)
                    node = parentNode.Children(i);
                    return
                end
            end

            node = uitreenode(parentNode, 'Text', text, 'UserData', 1);
        end

        function updateResult(obj, result, updateReactants)
            % Project one semantic result case into the Results tab
            %
            % Args:
            %     result (struct): Semantic result case
            %     updateReactants (logical): Flag to refresh the reactants table
            obj.updateResultTabs(result);
            obj.updateMixtureProperties(result);
            obj.updateTables(result, updateReactants);
            obj.updateEquivalenceRatio(result);
            obj.updateUpstreamTurbulence(result);
            obj.updateError(result);
        end

        function updateResultTabs(obj, result)
            % Update result-tab visibility from the selected result content
            %
            % Args:
            %     result (struct): Semantic result case
            if obj.hasOutputState(result, 'statistics')
                obj.showTab('TurbulenceStatistics');
                obj.setVisible('Panel_parameters_2', 'on');
            else
                obj.hideTab('TurbulencestatisticsTab');
                obj.setVisible('Panel_parameters_2', 'off');
            end
        end

        function updateMixtureProperties(obj, result)
            mixtureIndex = 0;

            for i = 1:numel(result.outputStates)
                outputState = result.outputStates(i);

                if ~strcmp(outputState.type, 'mixture') || ~isfield(result, outputState.field)
                    continue
                end

                mixture = result.(outputState.field);

                if ~obj.isMixture(mixture)
                    continue
                end

                mixtureIndex = mixtureIndex + 1;

                if mixtureIndex > 5
                    return
                end

                suffix = ['_', num2str(mixtureIndex)];
                obj.updateCommonProperties(mixture, suffix);
                obj.updateFlowProperties(mixture, suffix);
                obj.updateRocketProperties(mixture, suffix);
                obj.updateObliqueProperties(mixture, suffix);
            end
        end

        function updateCommonProperties(obj, mixture, suffix)
            propertyMap = { ...
                'text_T', @(mix) temperature(mix); ...
                'text_p', @(mix) pressure(mix); ...
                'text_r', @(mix) density(mix); ...
                'text_h', @(mix) enthalpy_mass(mix); ...
                'text_e', @(mix) intEnergy_mass(mix); ...
                'text_cp', @(mix) cp_mass(mix); ...
                'text_s', @(mix) entropy_mass(mix); ...
                'text_gamma', @(mix) adiabaticIndex(mix); ...
                'text_W', @(mix) MolecularWeight(mix); ...
                'text_sound', @(mix) soundspeed(mix)};

            for i = 1:size(propertyMap, 1)
                obj.setNumericValue([propertyMap{i, 1}, suffix], ...
                    obj.evaluateProperty(propertyMap{i, 2}, mixture, Inf));
            end
        end

        function updateFlowProperties(obj, mixture, suffix)
            if ~isprop(mixture, 'uShock')
                return
            end

            if strcmp(suffix, '_1')
                velocity = obj.evaluateProperty(@(mix) velocity_relative(mix), mixture, 0);
            else
                velocity = obj.scalarProperty(mixture, 'uShock', 0);
            end

            obj.setNumericValue(['text_u', suffix], velocity);
            obj.setNumericValue(['text_M', suffix], velocity / soundspeed(mixture));
        end

        function updateRocketProperties(obj, mixture, suffix)
            obj.setNumericValue(['text_Aratio', suffix], obj.scalarProperty(mixture, 'areaRatio', 0));
            obj.setNumericValue(['text_Cstar', suffix], obj.scalarProperty(mixture, 'cstar', 0));
            obj.setNumericValue(['text_Ivac', suffix], obj.scalarProperty(mixture, 'I_vac', 0));
            obj.setNumericValue(['text_Isp', suffix], obj.scalarProperty(mixture, 'I_sp', 0));
        end

        function updateObliqueProperties(obj, mixture, suffix)
            obj.setNumericValue(['text_beta_min', suffix], obj.scalarProperty(mixture, 'betaMin', 0));
            obj.setNumericValue(['text_beta', suffix], obj.scalarProperty(mixture, 'beta', 0));
            obj.setNumericValue(['text_theta', suffix], obj.scalarProperty(mixture, 'theta', 0));
        end

        function updateTables(obj, result, updateReactants)
            initialMixture = obj.initialMixture(result);
            selectedMixture = obj.selectedMixture(result);

            if updateReactants && obj.isMixture(initialMixture)
                obj.updateReactantsTable(result);
            end

            if obj.isMixture(selectedMixture)
                obj.updateProductsTable(selectedMixture);
            end
        end

        function updateReactantsTable(obj, result)
            if ~isfield(result, 'UITable_R_Data') || isempty(result.UITable_R_Data)
                return
            end

            data = result.UITable_R_Data;
            obj.app.UITable_R2.Data = data(:, 1:3);
        end

        function updateProductsTable(obj, mixture)
            species = obj.mixtureSpecies(mixture);
            [species, moles, fractions] = obj.sortedComposition(mixture, species);
            obj.app.UITable_P.Data = table2cell(table(species, moles, fractions));
        end

        function [species, moles, fractions, varargout] = sortedComposition(obj, mixture, species, varargin)
            species = species(:);
            compositionSpecies = obj.mixtureSpecies(mixture);
            [~, index] = ismember(species, compositionSpecies);
            valid = index > 0;
            species = species(valid);
            index = index(valid);

            for i = 1:numel(varargin)
                varargin{i} = varargin{i}(valid);
            end

            fractions = mixture.Xi(index);
            fractions = fractions(:);
            moles = fractions .* mixture.N;
            [fractions, order] = sort(fractions, 'descend');
            moles = moles(order);
            species = species(order);

            varargout = cell(1, numel(varargin));

            for i = 1:numel(varargin)
                varargout{i} = varargin{i}(order);
            end
        end

        function updateEquivalenceRatio(obj, result)
            mixture = obj.initialMixture(result);

            if ~obj.isMixture(mixture) || isempty(mixture.equivalenceRatio)
                obj.setTextValue('edit_phi2', '-');
                obj.setTextValue('edit_phi3', '-');
                return
            end

            equivalenceRatio = sprintf('%.5g', round(mixture.equivalenceRatio, 5));
            obj.setTextValue('edit_phi2', equivalenceRatio);
            obj.setTextValue('edit_phi3', equivalenceRatio);
        end

        function updateUpstreamTurbulence(obj, result)
            mixture = obj.initialMixture(result);
            selectedMixture = obj.selectedMixture(result);

            if ~obj.isMixture(selectedMixture) || isempty(selectedMixture.lia)
                return
            end

            obj.setTextValue('edit_eta', sprintf('%.5g', round(mixture.eta, 5)));
            obj.setTextValue('edit_chi', sprintf('%.5g', round(mixture.chi, 5)));
            obj.setTextValue('edit_etaVorticity', sprintf('%.5g', round(mixture.etaVorticity, 5)));
            obj.updateTurbulenceStatistics(result);
        end

        function updateTurbulenceStatistics(obj, result)
            statistics = obj.turbulenceStatistics(result);

            if isempty(statistics)
                return
            end

            propertyMap = { ...
                'text_K_2', 'K'; ...
                'text_R11_2', 'R11'; ...
                'text_RTT_2', 'RTT'; ...
                'text_Enstrophy_2', 'enstrophy'; ...
                'text_EnstrophyTT_2', 'enstrophyTT'; ...
                'text_Ka_2', 'Ka'; ...
                'text_R11a_2', 'R11a'; ...
                'text_RTTa_2', 'RTTa'; ...
                'text_Kr_2', 'Kr'; ...
                'text_R11r_2', 'R11r'; ...
                'text_RTTr_2', 'RTTr'; ...
                'text_Kolmogorov_2', 'kolmogorovLengthRatio'};

            for i = 1:size(propertyMap, 1)
                obj.setNumericValue(propertyMap{i, 1}, ...
                    obj.statisticsValue(statistics, propertyMap{i, 2}, 0));
            end
        end

        function updateError(obj, result)
            errors = [];

            for i = 1:numel(result.outputStates)
                outputState = result.outputStates(i);

                if ~strcmp(outputState.type, 'mixture') || ~isfield(result, outputState.field)
                    continue
                end

                mixture = result.(outputState.field);

                if ~obj.isMixture(mixture)
                    continue
                end

                if strcmpi(result.ProblemType, 'TP') && isprop(mixture, 'errorMoles') && ~isempty(mixture.errorMoles)
                    errors(end + 1) = mixture.errorMoles; %#ok<AGROW>
                elseif isprop(mixture, 'errorProblem') && ~isempty(mixture.errorProblem)
                    errors(end + 1) = mixture.errorProblem; %#ok<AGROW>
                end
            end

            if isempty(errors)
                obj.setNumericValue('text_error_problem', 0);
                return
            end

            obj.setNumericValue('text_error_problem', max(errors));
        end

        function writeSelectedCaseReport(obj, result, caseLabel)
            if isempty(obj.app)
                return
            end

            lines = obj.selectedCaseReport(result, caseLabel);
            message = strjoin(lines, newline);

            if isprop(obj.app, 'consolePanel') && ~isempty(obj.app.consolePanel)
                obj.app.consolePanel.setOutput(message);
            elseif isprop(obj.app, 'Console_text') && isprop(obj.app.Console_text, 'Value')
                obj.app.Console_text.Value = message;
            end
        end

        function lines = selectedCaseReport(obj, result, caseLabel)
            setup = obj.setupFromResult(result);
            summaryRows = { ...
                'Case', obj.caseText(caseLabel); ...
                'Problem type', obj.textValue(result.ProblemType, '-')};

            if isempty(fieldnames(setup))
                summaryRows(end + 1, :) = {'Setup', 'not available'};
                lines = obj.keyValueLines(summaryRows);
                return
            end

            if obj.hasEquivalenceRatio(result, setup)
                summaryRows(end + 1, :) = {'Equivalence ratio', obj.equivalenceRatioText(result, setup)};
            end

            lines = obj.keyValueLines(summaryRows);
            lines = [lines; obj.statePropertiesLines( ...
                'Reactant state', setup.reactantStateProperties, result)];
            lines = [lines; obj.productStateLines(setup, result)];
            lines = [lines; obj.inputPairLines( ...
                'Additional reactant inputs', setup.additionalInputsReactants, result)];
            lines = [lines; obj.inputPairLines( ...
                'Additional product inputs', setup.additionalInputsProducts, result)];
            lines = [lines; obj.flagsLines(setup.flags)];
            lines = [lines; obj.reactantSpeciesLines(result, setup)];
            lines = [lines; obj.productSpeciesLines(result, setup)];
        end

        function setup = setupFromResult(~, result)
            setup = struct();

            if isfield(result, 'setup') && isstruct(result.setup)
                setup = result.setup;
            end
        end

        function lines = statePropertiesLines(obj, title, stateProperties, result)
            lines = {};

            if isempty(stateProperties)
                return
            end

            rows = cell(numel(stateProperties), 2);

            for i = 1:numel(stateProperties)
                value = obj.caseValue(stateProperties(i).value, result);
                rows(i, :) = { ...
                    obj.statePropertyLabel(stateProperties(i).name), ...
                    obj.formatValue(value)};
            end

            lines = obj.tableLines(title, {'Property', 'Value'}, rows);
        end

        function lines = productStateLines(obj, setup, result)
            lines = {};

            if ~obj.hasProductStateConstraint(setup)
                return
            end

            problemType = upper(char(setup.requestedProblemType));

            switch problemType
                case 'HP'
                    lines = obj.mirroredProductStateLines('enthalpy', 'pressure', setup, result);
                case 'SP'
                    lines = obj.mirroredProductStateLines('entropy', 'pressure', setup, result);
                case 'EV'
                    lines = obj.mirroredProductStateLines('internal energy', 'volume', setup, result);
                case 'SV'
                    lines = obj.mirroredProductStateLines('entropy', 'volume', setup, result);
                otherwise
                    lines = obj.statePropertiesLines('Product state', setup.productStateProperties, result);
            end
        end

        function value = hasProductStateConstraint(~, setup)
            value = isfield(setup, 'constraint') && ~isempty(setup.constraint);
        end

        function lines = mirroredProductStateLines(obj, mirroredProperty, mechanicalProperty, setup, result)
            rows = {obj.statePropertyLabel(mirroredProperty), 'reactants'};
            value = obj.statePropertyValue(setup.productStateProperties, mechanicalProperty, result);

            if ~isempty(value)
                rows(end + 1, :) = {obj.statePropertyLabel(mechanicalProperty), ...
                    obj.formatValue(value)};
            end

            lines = obj.tableLines('Product state', {'Property', 'Value'}, rows);
        end

        function value = statePropertyValue(obj, stateProperties, name, result)
            value = [];

            if isempty(stateProperties)
                return
            end

            index = find(strcmpi({stateProperties.name}, name), 1);

            if isempty(index)
                return
            end

            value = obj.caseValue(stateProperties(index).value, result);
        end

        function lines = inputPairLines(obj, title, inputs, result)
            lines = {};

            if isempty(inputs)
                return
            end

            rows = {};

            for i = 1:2:numel(inputs)
                if i + 1 > numel(inputs)
                    continue
                end

                if strcmpi(inputs{i}, 'equivalenceRatio')
                    continue
                end

                value = obj.caseValue(inputs{i + 1}, result);
                rows(end + 1, :) = {obj.textValue(inputs{i}, '-'), ...
                    obj.formatValue(value)}; %#ok<AGROW>
            end

            lines = obj.tableLines(title, {'Input', 'Value'}, rows);
        end

        function lines = flagsLines(obj, flags)
            lines = {};

            if isempty(flags)
                return
            end

            names = fieldnames(flags);

            if isempty(names)
                return
            end

            lines = {'Options:'};

            for i = 1:numel(names)
                line = obj.flagReportLine(flags, names{i});

                if isempty(line)
                    continue
                end

                lines{end + 1} = ['  ', line]; %#ok<AGROW>
            end

            if isscalar(lines)
                lines = {};
                return
            end

            lines = lines(:);
        end

        function lines = reactantSpeciesLines(obj, result, setup)
            lines = {};
            mixture = obj.initialMixture(result);

            if obj.isMixture(mixture)
                species = mixture.listSpecies(:);
                amount = mixture.quantity(:);
                rows = cell(numel(species), 4);

                for i = 1:numel(species)
                    role = obj.speciesRole(mixture, species{i});
                    temperature = obj.reactantSpeciesTemperature(result, setup, species{i});
                    rows(i, :) = {species{i}, obj.formatValue(amount(i)), ...
                        role, obj.formatValue(temperature)};
                end

                lines = obj.tableLines(sprintf('Reactants (%d)', numel(species)), ...
                    {'Species', 'Moles', 'Type', 'T [K]'}, rows);
                return
            end

            if ~isfield(setup, 'reactants') || ~isfield(setup.reactants, 'species') ...
                    || isempty(setup.reactants.species)
                return
            end

            rows = cell(numel(setup.reactants.species), 4);

            for i = 1:numel(setup.reactants.species)
                rows(i, :) = {setup.reactants.species{i}, ...
                    obj.formatValue(setup.reactants.amount(i)), '-', '-'};
            end

            lines = obj.tableLines(sprintf('Reactants (%d)', numel(setup.reactants.species)), ...
                {'Species', 'Moles', 'Type', 'T [K]'}, rows);
        end

        function lines = productSpeciesLines(obj, result, setup)
            species = {};

            if isfield(result, 'listProducts') && ~isempty(result.listProducts)
                species = result.listProducts;
            elseif isfield(result, 'listSpecies') && ~isempty(result.listSpecies)
                species = result.listSpecies;
            elseif isfield(setup, 'products') && isfield(setup.products, 'species')
                species = setup.products.species;
            end

            if isempty(species)
                lines = {};
                return
            end

            lines = { ...
                ''; ...
                sprintf('Products (%s, %d)', obj.productSetText(setup), numel(species)); ...
                ['  ', obj.listText(species, 18, false)]};
        end

        function value = caseText(obj, caseLabel)
            value = obj.textValue(caseLabel, 'Case');

            if startsWith(value, 'Case ', 'IgnoreCase', true)
                value = regexprep(value, '^Case\s+', '');
                return
            end
        end

        function lines = keyValueLines(obj, rows)
            rows = obj.tableText(rows);
            lines = cell(size(rows, 1), 1);

            for i = 1:size(rows, 1)
                lines{i} = [rows{i, 1}, ' ', rows{i, 2}];
            end
        end

        function lines = tableLines(obj, title, headers, rows)
            if isempty(rows)
                lines = {};
                return
            end

            headers = obj.tableText(headers);
            rows = obj.tableText(rows);
            widths = obj.tableColumnWidths(headers, rows);
            lines = cell(size(rows, 1) + 4, 1);
            lines{1} = '';
            lines{2} = title;
            lines{3} = obj.tableRow(headers, widths);
            lines{4} = obj.tableSeparator(widths);

            for i = 1:size(rows, 1)
                lines{i + 4} = obj.tableRow(rows(i, :), widths);
            end
        end

        function values = tableText(obj, values)
            if isstring(values)
                values = cellstr(values);
            elseif ~iscell(values)
                values = {values};
            end

            for i = 1:numel(values)
                values{i} = obj.textValue(values{i}, '-');
            end
        end

        function widths = tableColumnWidths(~, headers, rows)
            widths = zeros(1, numel(headers));

            for i = 1:numel(headers)
                widths(i) = numel(headers{i});
            end

            for i = 1:size(rows, 1)
                for j = 1:size(rows, 2)
                    widths(j) = max(widths(j), numel(rows{i, j}));
                end
            end
        end

        function line = tableRow(~, values, widths)
            parts = cell(1, numel(values));

            for i = 1:numel(values)
                parts{i} = sprintf(['%-', num2str(widths(i)), 's'], values{i});
            end

            line = ['  ', strjoin(parts, '  ')];
        end

        function line = tableSeparator(~, widths)
            parts = cell(1, numel(widths));

            for i = 1:numel(widths)
                parts{i} = repmat('-', 1, widths(i));
            end

            line = ['  ', strjoin(parts, '  ')];
        end

        function value = caseValue(obj, value, result)
            if ~(isnumeric(value) || islogical(value)) || numel(value) <= 1
                return
            end

            index = obj.caseIndex(result);
            index = min(index, numel(value));
            value = value(index);
        end

        function index = caseIndex(~, result)
            index = 1;

            if isfield(result, 'caseIndex') && ~isempty(result.caseIndex)
                index = result.caseIndex(1);
            end

            index = max(1, round(index));
        end

        function value = productSetText(obj, setup)
            value = 'All';

            if isfield(setup, 'flags') && isfield(setup.flags, 'frozenChemistry') ...
                    && setup.flags.frozenChemistry
                value = 'Frozen';
                return
            end

            if ~isfield(setup, 'products') || ~isfield(setup.products, 'selection')
                return
            end

            value = obj.textValue(setup.products.selection, 'All');

            if strcmpi(value, 'Default')
                value = 'All';
            end
        end

        function label = statePropertyLabel(~, name)
            switch lower(char(name))
                case 'temperature'
                    label = 'temperature [K]';
                case 'pressure'
                    label = 'pressure [bar]';
                case 'volume'
                    label = 'specific volume [m3/kg]';
                case {'entropy', 'entropyspecific'}
                    label = 'entropy [kJ/(kg K)]';
                case {'enthalpy', 'enthalpyspecific'}
                    label = 'enthalpy [kJ/kg]';
                case {'internal energy', 'internalenergy', 'internalenergyspecific'}
                    label = 'internal energy [kJ/kg]';
                otherwise
                    label = char(name);
            end
        end

        function value = hasEquivalenceRatio(obj, result, setup)
            mixture = obj.initialMixture(result);
            value = obj.isMixture(mixture) && ~isempty(mixture.equivalenceRatio);

            if value || ~isfield(setup, 'flags') || ~isfield(setup.flags, 'hasEquivalenceRatio')
                return
            end

            value = setup.flags.hasEquivalenceRatio;
        end

        function value = equivalenceRatioText(obj, result, setup)
            mixture = obj.initialMixture(result);

            if obj.isMixture(mixture) && ~isempty(mixture.equivalenceRatio)
                value = obj.formatValue(mixture.equivalenceRatio);
                return
            end

            if isfield(setup, 'additionalInputsReactants')
                inputs = setup.additionalInputsReactants;

                for i = 1:2:numel(inputs)
                    if strcmpi(inputs{i}, 'equivalenceRatio') && i + 1 <= numel(inputs)
                        value = obj.formatValue(obj.caseValue(inputs{i + 1}, result));
                        return
                    end
                end
            end

            value = '-';
        end

        function line = flagReportLine(obj, flags, name)
            line = '';
            value = flags.(name);

            if isempty(value) || strcmpi(name, 'hasEquivalenceRatio')
                return
            end

            if strcmpi(name, 'infiniteAreaChamber')

                if ~value
                    line = 'finite area chamber';
                end

                return
            end

            if strcmpi(name, 'rocketSubsonic')

                if ~isfield(flags, 'areaRatio') || ~flags.areaRatio
                    return
                end

                if value
                    line = 'subsonic area ratio';
                else
                    line = 'supersonic area ratio';
                end

                return
            end

            if islogical(value) && isscalar(value)

                if value
                    line = obj.flagLabel(name);
                end

                return
            end

            line = [obj.flagLabel(name), ': ', obj.formatValue(value)];
        end

        function label = flagLabel(~, name)
            switch name
                case 'frozenChemistry'
                    label = 'frozen chemistry';
                case 'ionizedSpecies'
                    label = 'ionized species';
                case 'idealAir'
                    label = 'ideal air';
                case 'printResults'
                    label = 'print results';
                case 'beta'
                    label = 'beta branch';
                case 'theta'
                    label = 'theta branch';
                case 'areaRatio'
                    label = 'area ratio';
                otherwise
                    label = name;
            end
        end

        function role = speciesRole(~, mixture, species)
            role = 'Reactant';

            if any(strcmp(mixture.listSpeciesFuel, species))
                role = 'Fuel';
            elseif any(strcmp(mixture.listSpeciesOxidizer, species))
                role = 'Oxidizer';
            elseif any(strcmp(mixture.listSpeciesInert, species))
                role = 'Inert';
            end
        end

        function temperature = reactantSpeciesTemperature(obj, result, setup, species)
            temperature = [];

            if isfield(result, 'UITable_R_Data') && size(result.UITable_R_Data, 2) >= 5
                data = result.UITable_R_Data;
                index = find(strcmp(data(:, 1), species), 1);

                if ~isempty(index)
                    if iscell(data)
                        temperature = data{index, 5};
                    else
                        temperature = data(index, 5);
                    end

                    return
                end
            end

            if isfield(setup, 'reactants') && isfield(setup.reactants, 'temperature') ...
                    && ~isempty(setup.reactants.temperature)
                index = find(strcmp(setup.reactants.species, species), 1);

                if ~isempty(index)
                    temperature = setup.reactants.temperature(index);
                    return
                end
            end

            mixture = obj.initialMixture(result);

            if obj.isMixture(mixture)
                temperature = mixture.T;
            end
        end

        function value = textValue(~, value, defaultValue)
            if nargin < 3
                defaultValue = '-';
            end

            if isempty(value)
                value = defaultValue;
            elseif ischar(value)
                value = char(value);
            elseif isstring(value)
                value = char(value);
            elseif isnumeric(value) || islogical(value)
                value = num2str(value);
            else
                value = defaultValue;
            end

            if isempty(strtrim(value))
                value = defaultValue;
            end
        end

        function text = formatValue(obj, value)
            if isempty(value)
                text = '-';
            elseif ischar(value)
                text = char(value);
            elseif isstring(value)
                if isscalar(value)
                    text = char(value);
                else
                    text = obj.listText(cellstr(value), 8);
                end
            elseif islogical(value)
                if isscalar(value)
                    text = obj.logicalText(value);
                else
                    text = obj.listText(num2cell(value), 8);
                end
            elseif isnumeric(value)
                text = obj.numericText(value);
            elseif iscell(value)
                text = obj.listText(value, 8);
            else
                text = class(value);
            end
        end

        function text = numericText(~, value)
            if isempty(value)
                text = '-';
                return
            end

            value = value(:)';
            maxItems = 8;
            count = numel(value);
            limit = min(count, maxItems);
            items = arrayfun(@(item) sprintf('%.6g', item), value(1:limit), ...
                'UniformOutput', false);

            if count == 1
                text = items{1};
                return
            end

            text = ['[', strjoin(items, ', ')];

            if count > limit
                text = [text, sprintf(', ... (%d total)', count)];
            end

            text = [text, ']'];
        end

        function text = listText(obj, value, maxItems, varargin)
            includeTotal = true;

            if nargin > 3
                includeTotal = varargin{1};
            end

            if ischar(value) || (isstring(value) && isscalar(value))
                text = char(value);
                return
            end

            if isstring(value)
                value = cellstr(value);
            end

            if ~iscell(value)
                text = obj.formatValue(value);
                return
            end

            value = value(:)';
            count = numel(value);
            limit = min(count, maxItems);
            items = cell(1, limit);

            for i = 1:limit
                items{i} = obj.formatValue(value{i});
            end

            text = strjoin(items, ', ');

            if includeTotal && count > limit
                text = [text, sprintf(', ... (%d total)', count)];
            elseif count > limit
                text = [text, ', ...'];
            end
        end

        function text = logicalText(~, value)
            if value
                text = 'true';
            else
                text = 'false';
            end
        end

        function updateCustomFigures(obj, results)
            delete(obj.app.Mixtures.Children);
            delete(obj.app.Variable_x.Children);
            delete(obj.app.Variable_y.Children);

            outputStates = results(1).outputStates;
            mixtureFields = {};

            for i = 1:numel(outputStates)
                if strcmp(outputStates(i).type, 'mixture')
                    mixtureFields{end + 1} = outputStates(i).field; %#ok<AGROW>
                end
            end

            obj.addTreeNodes(obj.app.Mixtures, mixtureFields);
            properties = obj.numericMixtureProperties(results(1), mixtureFields);
            obj.addTreeNodes(obj.app.Variable_x, properties);
            obj.addTreeNodes(obj.app.Variable_y, properties);
        end

        function properties = numericMixtureProperties(~, result, mixtureFields)
            properties = {};

            if isempty(mixtureFields)
                return
            end

            mixture = result.(mixtureFields{end});
            names = sort(fieldnames(mixture));

            for i = 1:numel(names)
                value = mixture.(names{i});

                if isfloat(value) && ~isempty(value)
                    properties{end + 1} = names{i}; %#ok<AGROW>
                end
            end
        end

        function addTreeNodes(~, parentNode, names)
            for i = 1:numel(names)
                uitreenode(parentNode, 'Text', names{i});
            end
        end

        function mixArrays = mixtureOutputArrays(~, results)
            outputStates = results(1).outputStates;
            mixArrays = {};

            for i = 1:numel(outputStates)
                outputState = outputStates(i);

                if ~strcmp(outputState.type, 'mixture')
                    continue
                end

                values = {};

                for j = 1:numel(results)
                    if isfield(results(j), outputState.field)
                        values{end + 1} = results(j).(outputState.field); %#ok<AGROW>
                    end
                end

                if ~isempty(values)
                    mixArrays{end + 1} = [values{:}]; %#ok<AGROW>
                end
            end
        end

        function names = checkedNodeText(obj, componentName, varargin)
            names = {};

            if ~isprop(obj.app, componentName)
                return
            end

            nodes = obj.app.(componentName).CheckedNodes;

            for i = 1:numel(nodes)
                text = nodes(i).Text;

                if nargin > 2 && strcmp(text, varargin{1})
                    continue
                end

                names{end + 1} = text; %#ok<AGROW>
            end
        end

        function config = currentPlotConfig(obj, xField, yFields)
            if isprop(obj.app, 'plotConfig') && ~isempty(obj.app.plotConfig)
                config = obj.app.plotConfig;
            else
                config = combustiontoolbox.utils.display.PlotConfig();
            end

            config.labelx = combustiontoolbox.utils.display.interpreterLabel(xField, config.label_type);
            config.labely = combustiontoolbox.utils.display.interpreterLabel(yFields, config.label_type);
        end

        function value = defaultPlotSettingsEnabled(obj)
            value = isprop(obj.app, 'DefaultsettingsCheckBox') && obj.app.DefaultsettingsCheckBox.Value;
        end

        function mixtures = mixtureSeries(~, results, field)
            mixtures = {};

            for i = 1:numel(results)
                if isfield(results(i), field)
                    mixtures{end + 1} = results(i).(field); %#ok<AGROW>
                end
            end
        end

        function applyLegend(obj, legendNames, numMixtures, numProperties, config)
            legendNames = legendNames(~cellfun(@isempty, legendNames));

            if numProperties + numMixtures > 2
                obj.app.UIAxes.YLabel.String = 'Multiple variables';
                combustiontoolbox.utils.display.setLegends(obj.app.UIAxes, flip(legendNames), 'config', config);
                obj.app.UIAxes.Legend.Visible = 'on';
                return
            end

            if ~isempty(obj.app.UIAxes.Legend)
                obj.app.UIAxes.Legend.Visible = 'off';
            end
        end

        function value = convertPlotUnits(~, fieldName, value, basis)
            switch lower(fieldName)
                case {'cp', 'cv', 'hf', 'ef', 'h', 'e', 'g', 's', 'det', 'dht', 'ds', 's0'}
                    value = value ./ basis;
                case {'mw', 'w'}
                    value = value * 1e3;
            end
        end

        function value = isCompositionProperty(~, fieldName)
            value = any(strcmpi(fieldName, {'Xi', 'Yi', 'Ni', 'moles', 'molarFraction', 'massFraction'}));
        end

        function mixture = initialMixture(obj, result)
            mixture = obj.mixtureByRole(result, 'initial');
        end

        function mixture = selectedMixture(~, result)
            mixture = [];

            for i = 1:numel(result.outputStates)
                outputState = result.outputStates(i);

                if strcmp(outputState.type, 'mixture') && ~strcmp(outputState.role, 'initial')
                    mixture = result.(outputState.field);
                    return
                end
            end
        end

        function mixture = mixtureByRole(~, result, role)
            mixture = [];

            for i = 1:numel(result.outputStates)
                outputState = result.outputStates(i);

                if strcmp(outputState.type, 'mixture') && strcmp(outputState.role, role)
                    mixture = result.(outputState.field);
                    return
                end
            end
        end

        function species = mixtureSpecies(~, mixture)
            species = mixture.chemicalSystem.listSpecies(:);

            if isempty(species)
                species = mixture.listSpecies(:);
            end
        end

        function value = evaluateProperty(~, evaluator, mixture, defaultValue)
            try
                value = evaluator(mixture);
            catch
                value = defaultValue;
            end
        end

        function value = scalarProperty(~, object, name, defaultValue)
            value = defaultValue;

            if isprop(object, name) && ~isempty(object.(name))
                value = object.(name);
            end
        end

        function statistics = turbulenceStatistics(~, result)
            statistics = [];

            if isfield(result, 'states') && isfield(result.states, 'statistics')
                statistics = result.states.statistics;
            elseif isfield(result, 'resultsLIA')
                statistics = result.resultsLIA;
            end
        end

        function value = statisticsValue(~, statistics, name, defaultValue)
            value = defaultValue;

            if isstruct(statistics) && isfield(statistics, name) && ~isempty(statistics.(name))
                value = statistics.(name);
            end
        end

        function value = hasOutputState(~, result, stateId)
            value = false;

            if ~isstruct(result) || ~isfield(result, 'outputStates')
                return
            end

            outputStates = result.outputStates;

            for i = 1:numel(outputStates)
                if strcmp(outputStates(i).id, stateId)
                    value = true;
                    return
                end
            end
        end

        function showTab(obj, tabTitle, tabProperty)
            if ~isempty(obj.app) && ismethod(obj.app, 'showTab')
                obj.app.showTab(tabTitle);
            end

            if nargin > 2
                obj.selectTab(tabProperty);
            end
        end

        function hideTab(obj, tabProperty)
            if ~isempty(obj.app) && isprop(obj.app, tabProperty) ...
                    && ismethod(obj.app, 'hideTab')
                obj.app.hideTab(obj.app.(tabProperty));
            end
        end

        function selectTab(obj, tabProperty)
            if isempty(obj.app) || ~isprop(obj.app, tabProperty) ...
                    || isempty(obj.app.(tabProperty)) ...
                    || ~isvalid(obj.app.(tabProperty)) ...
                    || isempty(obj.app.(tabProperty).Parent)
                return
            end

            obj.app.(tabProperty).Parent.SelectedTab = obj.app.(tabProperty);
        end

        function setNumericValue(obj, componentName, value)
            % Set a numeric component value when the component exists
            %
            % Args:
            %     componentName (char): App component name
            %     value: Value assigned to the component
            if isprop(obj.app, componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.numericValue(value);
                obj.app.(componentName).Value = value;
            end
        end

        function setTextValue(obj, componentName, value)
            % Set a text component value when the component exists
            %
            % Args:
            %     componentName (char): App component name
            %     value: Value assigned to the component
            if isprop(obj.app, componentName) && isprop(obj.app.(componentName), 'Value')
                obj.app.(componentName).Value = value;
            end
        end

        function setVisible(obj, componentName, value)
            % Set component visibility when the component exists
            %
            % Args:
            %     componentName (char): App component name
            %     value (char | logical): Visibility value
            if isprop(obj.app, componentName) && isprop(obj.app.(componentName), 'Visible')
                obj.app.(componentName).Visible = value;
            end
        end

        function value = isMixture(~, value)
            value = isa(value, 'combustiontoolbox.core.Mixture') && ~isempty(value);
        end

        function value = numericValue(~, value)
            if isempty(value)
                value = 0;
                return
            end

            if ischar(value) || isstring(value)
                value = str2double(value);
            end

            if ~isnumeric(value) || isempty(value)
                value = 0;
                return
            end

            value = double(value(1));

            if ~isreal(value) || ~isfinite(value)
                value = 0;
            end
        end

    end
end

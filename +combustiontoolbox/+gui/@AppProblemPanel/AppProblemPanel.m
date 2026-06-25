classdef AppProblemPanel < handle
    % Applies selected-problem metadata to GUI controls and feature panels.
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %     catalog (AppProblemCatalog): Problem definitions and metadata
    %     inputsPanel (AppProblemInputsPanel): Coupled input control updater
    %
    % Examples:
    %     * panel = combustiontoolbox.gui.AppProblemPanel(app, catalog);
    %     * panel.applySelectedProblem();

    properties (Access = private)
        app
        catalog
        inputsPanel
    end

    methods
        function obj = AppProblemPanel(app, catalog, inputsPanel)
            % AppProblemPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %     catalog (AppProblemCatalog): Problem definitions and metadata
            %     inputsPanel (AppProblemInputsPanel): Coupled input updater
            %
            % Returns:
            %     obj (AppProblemPanel): Initialized problem panel object
            if nargin < 1
                app = [];
            end

            if nargin < 2 || isempty(catalog)
                catalog = combustiontoolbox.gui.AppProblemCatalog();
            end

            if nargin < 3 || isempty(inputsPanel)
                inputsPanel = combustiontoolbox.gui.AppProblemInputsPanel(app);
            end

            obj.app = app;
            obj.catalog = catalog;
            obj.inputsPanel = inputsPanel;
        end

        function definition = applySelectedProblem(obj)
            % Apply the currently selected problem metadata to the GUI
            %
            % Returns:
            %     definition (struct): Applied problem definition
            problemType = obj.componentValue('ProblemType', 'TP');
            definition = obj.applyProblem(problemType);
        end

        function definition = applyProblem(obj, problemType)
            % Apply explicit problem metadata to the GUI
            %
            % Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     definition (struct): Applied problem definition
            problemType = obj.normalizeProblemType(problemType);
            definition = obj.catalog.get(problemType);

            obj.clearResults();

            if strcmp(definition.status, 'notIncluded')
                obj.showNotIncludedAlert();
                return
            end

            obj.applyProblemDefinition(definition);
            obj.updateReactantsTemperature();
        end

        function onRocketModelChanged(obj)
            % Apply rocket model controls
            obj.applyRocketModelControls();
        end

    end

    methods (Access = private)
        function problemType = normalizeProblemType(obj, problemType) %#ok<INUSD>
            problemType = char(problemType);

            if any(strcmp(problemType, {'ROCKET_IAC', 'ROCKET_FAC'}))
                problemType = 'ROCKET';
            end
        end

        function applyProblemDefinition(obj, definition)
            visibility = definition.visibility;
            labels = definition.labels;
            defaults = definition.defaultInputs;

            obj.applyIACFlag(visibility);
            obj.resetFeatureVisibility();
            obj.applyInputVisibility(visibility);
            obj.applyLabels(labels, visibility);
            obj.applyActiveFeatureVisibility(definition);
            obj.applyDefaultInputs(defaults);

            if obj.requiresMachUpdate(definition)
                obj.inputsPanel.updateMachOrVelocity('Mach');
            end
        end

        function applyIACFlag(obj, visibility)
            if visibility.FLAG_IAC
                obj.setVisible('FLAG_IAC', true);
            else
                obj.setVisible('FLAG_IAC', false);
                obj.setValue('FLAG_IAC', true);
            end

            obj.applyRocketModelControls();
        end

        function resetFeatureVisibility(obj)
            obj.setShockPanelVisible(false);
            obj.setObliquePanelVisible(false);
            obj.setRocketPanelVisible(false);
            obj.setShockTurbulencePanelVisible(false);
        end

        function applyInputVisibility(obj, visibility)
            obj.setVisible('PP1', visibility.PP1);
            obj.setVisible('PP2', visibility.PP2);
            obj.setVisible('PP3', visibility.PP3);
            obj.setVisible('PP4', visibility.PP4);
            obj.setVisible('PP5', visibility.PP5);
            obj.setVisible('PP6', visibility.PP6);
            obj.setVisible('PR3', visibility.PR3);
            obj.setVisible('PR4', visibility.PR4);
            obj.setVisible('PR5', visibility.PR5);
            obj.setVisible('text_P1', visibility.text_P1);
            obj.setVisible('AdditionalconstraintsPanel', visibility.AdditionalconstraintsPanel);
            obj.setVisible('text_RP', visibility.text_RP);
            obj.setVisible('text_R2', visibility.text_R2);
            obj.setVisible('text_P2', visibility.text_P2);
            obj.setVisible('text_RP4', visibility.text_RP4);
            obj.setVisible('text_RP5', visibility.text_RP5);
        end

        function applyLabels(obj, labels, visibility)
            obj.setText('text_P1', labels.text_P1);
            obj.setText('text_RP2', labels.text_RP2);
            obj.setText('text_RP', labels.text_RP);
            obj.setText('text_R2', labels.text_R2);
            obj.setText('text_P2', labels.text_P2);
            obj.setText('text_RP3', labels.text_RP3);
            obj.setText('text_RP4', labels.text_RP4);
            obj.setText('text_RP5', labels.text_RP5);
            obj.setText('text_RP1_2', labels.text_RP1_2);
            obj.setText('text_RP2_2', labels.text_RP2_2);
            obj.setText('text_RP3_2', labels.text_RP3_2);
            obj.setTitle('AdditionalconstraintsPanel', labels.additionalPanelTitle);

            if ~visibility.text_RP4
                obj.setVisible('text_RP4', false);
            end

            if ~visibility.text_RP5
                obj.setVisible('text_RP5', false);
            end
        end

        function applyActiveFeatureVisibility(obj, definition)
            visibility = definition.visibility;

            if visibility.shocks
                obj.setShockPanelVisible(true);
            end

            if visibility.oblique
                obj.setObliquePanelVisible(true);
            end

            if visibility.rocket
                obj.setRocketPanelVisible(true);
            end

            if visibility.shockTurbulence
                obj.setShockTurbulencePanelVisible(true, obj.shockTurbulenceModel(definition.id));
            end
        end

        function applyDefaultInputs(obj, defaults)
            fields = fieldnames(defaults);

            for i = 1:numel(fields)
                componentName = fields{i};
                value = defaults.(componentName);

                if ~obj.hasComponent(componentName)
                    continue
                end

                if ~isempty(value)
                    obj.setValue(componentName, value);
                elseif obj.componentIsVisible(componentName)
                    obj.setValue(componentName, '');
                end
            end
        end

        function value = requiresMachUpdate(obj, definition) %#ok<INUSD>
            value = strcmp(definition.labels.text_RP4, 'Mach number [-]') && ...
                ~isempty(definition.defaultInputs.PR4);
        end

        function updateReactantsTemperature(obj)
            if isempty(obj.componentData('UITable_R', {}))
                return
            end

            value = combustiontoolbox.gui.AppInput.parseValue(obj.componentValue('PR1', []), 'first');
            data = obj.componentData('UITable_R', {});
            data(:, 5) = repmat({value}, size(data(:, 5)));
            obj.setData('UITable_R', data);
        end

        function clearResults(obj)
            for i = 1:5
                suffix = sprintf('_%d', i);
                obj.updateCommonResults(0, suffix);
            end

            for i = 2:5
                suffix = sprintf('_%d', i);
                obj.updateRocketResults(0, suffix);
            end

            obj.setValue('text_error_problem', 0);
        end

        function updateCommonResults(obj, value, suffix)
            names = {'text_T', 'text_p', 'text_r', 'text_h', 'text_e', ...
                'text_cp', 'text_s', 'text_gamma', 'text_W', 'text_sound', ...
                'text_u', 'text_M'};

            for i = 1:numel(names)
                obj.setValue([names{i}, suffix], value);
            end
        end

        function updateRocketResults(obj, value, suffix)
            names = {'text_Aratio', 'text_Cstar', 'text_Ivac', 'text_Isp'};

            for i = 1:numel(names)
                obj.setValue([names{i}, suffix], value);
            end
        end

        function applyRocketModelControls(obj)
            if ~obj.componentValue('FLAG_IAC', true)
                obj.setText('text_P1', 'FAC');
                obj.setVisible('text_P1', true);
                obj.setVisible('text_RP1_2', true);
                obj.setVisible('text_RP2_2', false);
                obj.setVisible('PP1', true);
                obj.setVisible('PP2', false);
                obj.setValue('PP1', '');
                obj.setValue('PP2', '');
                obj.setVisible('Panel_extra_5', true);
                obj.setText('text_Products', 'Injector');
                obj.setText('text_Products_3', 'Outlet Chamber');
                obj.setText('text_Products_4', 'Throat');
                obj.setText('text_Products_5', 'Exit');
                return
            end

            obj.setText('text_P1', 'Products');
            obj.setVisible('text_P1', false);
            obj.setVisible('text_RP1_2', false);
            obj.setVisible('text_RP2_2', false);
            obj.setVisible('PP1', false);
            obj.setVisible('PP2', false);
            obj.setValue('PP1', '2500');
            obj.setValue('PP2', '1');
            obj.setVisible('Panel_extra_5', false);
            obj.setText('text_Products', 'Outlet Chamber');
            obj.setText('text_Products_3', 'Throat');
            obj.setText('text_Products_4', 'Exit');
            obj.setText('text_Products_5', 'Exit');
        end

        function setShockPanelVisible(obj, value)
            names = {'text_u', 'text_u_1', 'text_u_2', ...
                'text_M', 'text_M_1', 'text_M_2'};
            obj.setComponentsVisible(names, value);
        end

        function setObliquePanelVisible(obj, value)
            names = {'text_beta_min', 'text_beta_min_2', ...
                'text_beta', 'text_beta_2', 'text_theta', 'text_theta_2'};
            obj.setComponentsVisible(names, value);
        end

        function setRocketPanelVisible(obj, value)
            if value
                delta_pos_x = 77;
                panelPosition = obj.app.defaultLayout.Panel_parameters.Position;
                panelPosition(1) = panelPosition(1) - delta_pos_x;
                panelPosition(3) = panelPosition(3) + 3 * delta_pos_x;
                obj.setPosition('Panel_parameters', panelPosition);
                obj.applyRocketModelControls();
            else
                if isfield(obj.app.defaultLayout, 'Panel_parameters')
                    obj.setPosition('Panel_parameters', obj.app.defaultLayout.Panel_parameters.Position);
                end

                obj.setText('text_Products', 'Products');
                obj.setVisible('Panel_extra_5', false);
            end

            names = {'text_Aratio', 'text_Aratio_2', 'text_Cstar', ...
                'text_Ivac', 'text_Isp', 'Panel_extra_3', 'Panel_extra_4'};
            obj.setComponentsVisible(names, value);
        end

        function setShockTurbulencePanelVisible(obj, value, varargin)
            RP1_2_text = 'Area ratio A_c/A_t';
            RP2_2_text = 'Mass flux [kg/s]';
            RP3_2_text = 'etaVorticity [-]';
            P2_value = value;
            PP6_value = '';

            if value
                value = true;
                P1_text = 'Upstream turbulence';
                PP1_value = '';
                PP2_value = '';
                obj.setShockTurbulenceModel(varargin{:});
                obj.showTurbulenceStatisticsTab();
            else
                value = false;
                P1_text = 'Products';
                PP1_value = '2500';
                PP2_value = '1';
                obj.hideTurbulenceStatisticsTab();
            end

            if nargin == 3
                switch varargin{1}
                    case {'VORTICAL', 'ACOUSTIC'}
                        value = false;
                        P2_value = false;
                    case 'VORTICAL_ENTROPIC'
                        RP1_2_text = 'chi [-]';
                        PP1_value = '-0.1';
                        P2_value = false;
                    case 'COMPRESSIBLE'
                        RP1_2_text = 'eta [-]';
                        PP1_value = '0.1';
                        RP2_2_text = 'chi [-]';
                        PP2_value = '0';
                        RP3_2_text = 'etaVorticity [-]';
                        PP6_value = '0.04';
                end
            end

            obj.setText('text_P1', P1_text);
            obj.setVisible('text_P1', value);
            obj.setVisible('text_RP1_2', value);
            obj.setText('text_RP1_2', RP1_2_text);
            obj.setVisible('text_RP2_2', P2_value);
            obj.setText('text_RP2_2', RP2_2_text);
            obj.setVisible('text_RP3_2', P2_value);
            obj.setText('text_RP3_2', RP3_2_text);
            obj.setVisible('PP1', value);
            obj.setVisible('PP2', P2_value);
            obj.setVisible('PP6', P2_value);
            obj.setValue('PP1', PP1_value);
            obj.setValue('PP2', PP2_value);
            obj.setValue('PP6', PP6_value);
        end

        function setShockTurbulenceModel(obj, varargin)
            if isempty(varargin) || ~obj.hasComponent('shockTurbulenceSolver')
                return
            end

            obj.app.shockTurbulenceSolver.setShockTurbulenceModel(varargin{1});
        end

        function showTurbulenceStatisticsTab(obj)
            if ismethod(obj.app, 'showTab')
                obj.app.showTab('TurbulenceStatistics');
            end
        end

        function hideTurbulenceStatisticsTab(obj)
            if obj.hasComponent('TurbulencestatisticsTab') && ismethod(obj.app, 'hideTab')
                obj.app.hideTab(obj.app.TurbulencestatisticsTab);
            end
        end

        function model = shockTurbulenceModel(obj, problemType) %#ok<INUSD>
            model = strrep(problemType, 'SHOCKTURBULENCE_', '');
        end

        function showNotIncludedAlert(obj)
            if obj.hasComponent('UIFigure')
                uialert(obj.app.UIFigure, 'Problem not included yet. Sorry for the inconvenience.', 'Error')
            end
        end

        function value = componentValue(obj, componentName, defaultValue)
            value = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.app.(componentName).Value;
            end
        end

        function data = componentData(obj, componentName, defaultValue)
            data = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Data')
                data = obj.app.(componentName).Data;
            end
        end

        function setValue(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.assignmentValue(obj.app.(componentName), value);
                obj.app.(componentName).Value = value;
            end
        end

        function setData(obj, componentName, data)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Data')
                obj.app.(componentName).Data = data;
            end
        end

        function setText(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Text')
                obj.app.(componentName).Text = value;
            end
        end

        function setTitle(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Title')
                obj.app.(componentName).Title = value;
            end
        end

        function setPosition(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Position')
                obj.app.(componentName).Position = value;
            end
        end

        function setVisible(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Visible')
                obj.app.(componentName).Visible = obj.visibleState(value);
            end
        end

        function setComponentsVisible(obj, componentNames, value)
            for i = 1:numel(componentNames)
                obj.setVisible(componentNames{i}, value);
            end
        end

        function value = componentIsVisible(obj, componentName)
            value = false;

            if ~obj.hasComponent(componentName) || ~isprop(obj.app.(componentName), 'Visible')
                return
            end

            visible = obj.app.(componentName).Visible;

            if islogical(visible)
                value = visible;
            elseif isnumeric(visible)
                value = logical(visible);
            else
                value = strcmpi(char(visible), 'on');
            end
        end

        function value = hasComponent(obj, componentName)
            value = isobject(obj.app) && isprop(obj.app, componentName) && ~isempty(obj.app.(componentName));
        end

        function value = visibleState(obj, flag) %#ok<INUSD>
            if flag
                value = 'on';
            else
                value = 'off';
            end
        end

        function value = assignmentValue(~, component, value)
            if ~isnumeric(component.Value)
                return
            end

            if isempty(value)
                value = 0;
            elseif ischar(value) || isstring(value)
                numericValue = str2double(value);

                if isnan(numericValue)
                    value = 0;
                else
                    value = numericValue;
                end
            elseif ~isnumeric(value)
                value = 0;
            end

            if isnumeric(value) && ~isempty(value)
                value = value(1);
            end
        end

    end
end

classdef uipreferences < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        preferences_UIFigure     matlab.ui.Figure
        apply_button             matlab.ui.control.Button
        cancel_button            matlab.ui.control.Button
        ok_button                matlab.ui.control.Button
        tree                     matlab.ui.container.Tree
        combustion_toolbox_node  matlab.ui.container.TreeNode
        general_node             matlab.ui.container.TreeNode
        tuning_parameters_node   matlab.ui.container.TreeNode
        tuning_parameters_ct_equil_tptv_node     matlab.ui.container.TreeNode
        tuning_parameters_ct_equil_hpspevsv_node matlab.ui.container.TreeNode
        tuning_parameters_ct_sd_node             matlab.ui.container.TreeNode
        tuning_parameters_ct_rocket_node         matlab.ui.container.TreeNode
        miscellaneous_node   matlab.ui.container.TreeNode
        ContextMenu          matlab.ui.container.ContextMenu
        SnapshotMenu         matlab.ui.container.Menu
    end
    
    properties (Access = public)
        x0_panel_right = 206;
        y0_panel_right = 428;
        height_panel_0 = 38;
        delta_x = 9;
        delta_y = -12;
        width_amplification = 1.15;
        height_amplification = 1;
        width_box = 60;
        height_box = 22;
        width_right = 521
        background_color = [0.9098 0.9098 0.8902]
        caller_app % Handle to caller app
        dynamic_components
    end

    % app creation and deletion
    methods (Access = public)

        % Construct app
        function app = uipreferences(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                create_components(app)

                % Register the app with App Designer
                registerApp(app, app.preferences_UIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.preferences_UIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.preferences_UIFigure)
        end
    end
    
    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            % Get caller_app
            if nargin > 1
                app.caller_app = varargin{1};
            else
                fprintf('Running uipreferences in read mode.\n');
                app.caller_app = App();
            end
            % Create right panel based on first node uitree
            fill_combustion_toolbox_node(app);
            % Expand uitree
            expand(app.tree);
        end
        
        % Selection changed function: Tree
        function tree_selection_changed(app, event)
            node_tag = event.Source.SelectedNodes.Tag;

            % Delete Childrens of dynamic_components
            delete_objects(app, 'dynamic_components');
            % Fill GUI
            switch lower(node_tag)
                case 'main'
                    fill_combustion_toolbox_node(app);
                case 'general'
                    fill_general_node(app);
                case 'tuning'
                    fill_tuning_parameters(app);
                case 'tuning_ct_equil_tptv'
                    fill_tuning_parameters_ct_equil_tptv_node(app);
                case 'tuning_ct_equil_hpspevsv'
                    fill_tuning_parameters_ct_equil_hpspevsv_node(app);
                case 'tuning_ct_sd'
                    fill_tuning_parameters_ct_sd_node(app);
                case 'tuning_ct_rocket'
                    fill_tuning_parameters_ct_rocket_node(app);
                case 'miscellaneous'
                    fill_miscellaneous_node(app);
            end
        end
        
        % Set parameters to default
        function cancel_button_pushed(app, event)
            % Assign default values
            app.caller_app.C = Constants();
            app.caller_app.TN = TuningProperties();
            app.caller_app.Misc = Miscellaneous();
            % Delete app
            delete(app)
        end

        % Set value in caller_app
        function set_value(app, event, varargin)
            Tag = get_tag(app, event);
            app.caller_app.(Tag{1}).(Tag{2}) = event.Value;
        end

        % Set value in caller_app
        function set_value_function(app, event, varargin)
            Tag = get_tag(app, event);
            app.caller_app.(Tag{1}).(Tag{2}) = str2func(event.Value);
        end

        % Get value from caller_app
        function value = get_tag(app, event)
            value = split(event.Source.Tag, filesep);
        end
        
        % Delete all objects from object
        function app = delete_objects(app, parent)
            if ~isempty(app.(parent))
                childrens = fieldnames(app.(parent));
                for i = 1:length(childrens)
                    delete(app.(parent).(childrens{i}));
                end
                app.(parent) = [];
            end
        end
        
        % Create title panel right
        function create_title(app, title)
            % Create title1_label
            app.dynamic_components.title1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.title1_label.BackgroundColor = [0.8 0.8 0.8];
            app.dynamic_components.title1_label.FontWeight = 'bold';
            app.dynamic_components.title1_label.Position = [206 448 app.width_right 22];
            app.dynamic_components.title1_label.Text = ['  ', title];
        end

        % Create label panel's children
        function [width_label, height_label] = create_panel_label(app, panel_name, obj_num, title)

            % Definitions
            label_name = [panel_name, 'label', sprintf('%d', obj_num)];
            % Get width and height
            [width_label, height_label] = get_width_height_string(app, title);
            % Create uilabel
            app.dynamic_components.(label_name) = uilabel(app.dynamic_components.(panel_name));
            app.dynamic_components.(label_name).HorizontalAlignment = 'left';
            app.dynamic_components.(label_name).Position = [app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 width_label, height_label];
            app.dynamic_components.(label_name).Text = title;
        end

        % Create editfield panel's children
        function create_panel_edit_field(app, panel_name, obj_num, tag, title, varargin)
            % Create label panel's children
            width_label = create_panel_label(app, panel_name, obj_num, title);

            % Definitions
            edit_name = [panel_name, 'edit_field', sprintf('%d', obj_num)];
            % Default
            limits = [1 Inf];
            display_format = '%d';
            type = 'numeric';
            subtype = [];
            % Unpack
            for i = 1:2:nargin-5
                switch lower(varargin{i})
                    case 'limits'
                        limits = varargin{i+1};
                    case 'display_format'
                        display_format = varargin{i+1};
                    case 'type'
                        type = varargin{i+1};
                    case 'subtype'
                        subtype = varargin{i+1};
                end
            end

            app.dynamic_components.(edit_name) = uieditfield(app.dynamic_components.(panel_name), type);
            if strcmpi(type, 'numeric')
                app.dynamic_components.(edit_name).Limits = limits;
                app.dynamic_components.(edit_name).ValueDisplayFormat = display_format;
                app.dynamic_components.(edit_name).Value = app.caller_app.(tag{1}).(tag{2});
                app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value, true);
            else
                switch lower(subtype)
                    case 'function'
                        app.dynamic_components.(edit_name).Value = func2str(app.caller_app.(tag{1}).(tag{2}));
                        app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value_function, true);
                end
            end
            
            app.dynamic_components.(edit_name).Position = [width_label + 2 * app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 app.width_box, app.height_box];
            app.dynamic_components.(edit_name).Tag = fullfile(tag{1}, tag{2});
        end
        
        % Create dropdown panel's children
        function create_panel_dropdown(app, panel_name, obj_num, tag, title, items, varargin)
            % Create label panel's children
            width_label = create_panel_label(app, panel_name, obj_num, title);

            % Definitions
            dropdown_name = [panel_name, 'dropdown', sprintf('%d', obj_num)];
            % Default
            subtype = [];
            % Unpack
            for i = 1:2:nargin-6
                switch lower(varargin{i})
                    case 'subtype'
                        subtype = varargin{i+1};
                end
            end

            app.dynamic_components.(dropdown_name) = uidropdown(app.dynamic_components.(panel_name));
            app.dynamic_components.(dropdown_name).Items = items;

            switch lower(subtype)
                case 'function'
                    app.dynamic_components.(dropdown_name).Value = func2str(app.caller_app.(tag{1}).(tag{2}));
                    app.dynamic_components.(dropdown_name).ValueChangedFcn = createCallbackFcn(app, @set_value_function, true);
                otherwise
                    app.dynamic_components.(dropdown_name).Value = app.caller_app.(tag{1}).(tag{2});
                    app.dynamic_components.(dropdown_name).ValueChangedFcn = createCallbackFcn(app, @set_value, true);
            end
            
            app.dynamic_components.(dropdown_name).Position = [width_label + 2 * app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 1.2 * app.width_box, app.height_box];
            app.dynamic_components.(dropdown_name).Tag = fullfile(tag{1}, tag{2});
        end

        function panel_name = create_panel(app, position, title, varargin)
            % Create default panel
            
            % Definitions
            delta_x = app.delta_x;
            delta_y = app.delta_y;
            [width_label, height_label] = get_width_height_string(app, title);
            obj_num = 1;
            scrollable = 'on';
            % Unpack
            for i = 1:2:nargin-3
                switch lower(varargin{i})
                    case {'num', 'obj_num'}
                        obj_num = varargin{i+1};
                    case 'delta_x'
                        delta_x = varargin{i+1};
                    case 'delta_y'
                        delta_y = varargin{i+1};
                    case 'width_label'
                        width_label = varargin{i+1};
                    case 'height_label'
                        height_label = varargin{i+1};
                    case 'scrollable'
                        scrollable = varargin{i+1};
                end
            end
            
            panel_name = ['panel', sprintf('%d', obj_num)];
            label_name = ['label', sprintf('%d', obj_num)];

            % Create Panel
            app.dynamic_components.(panel_name) = uipanel(app.preferences_UIFigure);
            app.dynamic_components.(panel_name).BackgroundColor = app.background_color;
            app.dynamic_components.(panel_name).Position = position;
            app.dynamic_components.(panel_name).Scrollable = scrollable;
            
            % Create Label
            app.dynamic_components.(label_name) = uilabel(app.preferences_UIFigure);
            app.dynamic_components.(label_name).BackgroundColor = app.background_color;
            app.dynamic_components.(label_name).Position = ...
                [app.dynamic_components.(panel_name).Position(1) + delta_x,...
                 app.dynamic_components.(panel_name).Position(2) + app.dynamic_components.(panel_name).Position(4) + delta_y,...
                 width_label, height_label];
            app.dynamic_components.(label_name).Text = title;
        end

        function [width, height] = get_width_height_string(app, title)
            % Get width (pixels) of a string

            % Create temporal uicontrol with
            app.dynamic_components.temp = uicontrol(app.preferences_UIFigure, 'Style', 'text', 'String', title, 'Visible', 'off');
            temp_position = get(app.dynamic_components.temp, 'Extent');
            width = app.width_amplification * temp_position(3);
            height = app.height_amplification * temp_position(4);
            % Delete
            delete(app.dynamic_components.temp)
            app.dynamic_components = rmfield(app.dynamic_components, 'temp');
        end

        function ContextMenuOpening(app, event)
            gui_SnapshotMenuSelected(app.preferences_UIFigure);
        end

    end

    % Component initialization
    methods (Access = private)
        
        function create_components(app)
            % Create UIFigure and components
            
            % Create preferences_UIFigure and hide until all components are created
            app.preferences_UIFigure = uifigure('Visible', 'off');
            app.preferences_UIFigure.Color = app.background_color;
            app.preferences_UIFigure.Position = [500 300 740 480];
            app.preferences_UIFigure.Name = 'UIPreferences';
            app.preferences_UIFigure.Resize = 'off';
            app.preferences_UIFigure.Icon = 'logo_uipreferences.png';

            % Create tree
            create_tree(app);
            % Create combustion_toolbox_node
            create_combustion_toolbox_node(app);
            % Create general_node
            create_general_node(app);
            % Create tuning_parameters_node
            create_tuning_parameters_node(app);
            % Create miscellaneous_node
            create_miscellaneous_node(app);
            % Create ok_button
            create_ok_button(app);
            % Create cancel_button
            create_cancel_button(app);
            % Create apply_button
            create_apply_button(app);

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.preferences_UIFigure);
            app.ContextMenu.ContextMenuOpeningFcn = createCallbackFcn(app, @ContextMenuOpening, true);
            % Create SnapshotMenu
            app.SnapshotMenu = uimenu(app.ContextMenu);
            app.SnapshotMenu.Text = 'Snapshot';
            % Assign app.ContextMenu
            app.preferences_UIFigure.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.preferences_UIFigure.Visible = 'on';
        end
        
        function create_tree(app)
            % Create UITree
            app.tree = uitree(app.preferences_UIFigure);
            app.tree.Position = [12 40 178 430];
            app.tree.SelectionChangedFcn = createCallbackFcn(app, @tree_selection_changed, true);
        end
        
        function create_combustion_toolbox_node(app)
            % Create combustion_toolbox_node
            app.combustion_toolbox_node = uitreenode(app.tree);
            app.combustion_toolbox_node.Text = 'Combustion Toolbox';
            app.combustion_toolbox_node.Tag = 'main';
        end
        
        function create_general_node(app)
            % Create general_node
            app.general_node = uitreenode(app.combustion_toolbox_node);
            app.general_node.Text = 'General';
            app.general_node.Tag = 'general';
        end

        function create_tuning_parameters_node(app)
            % Create tuning_parameters_node
            app.tuning_parameters_node = uitreenode(app.combustion_toolbox_node);
            app.tuning_parameters_node.Text = 'Tuning parameters';
            app.tuning_parameters_node.Tag = 'tuning';
            % Create subnodes
            create_tuning_parameters_ct_equil_tptv_node(app);
            create_tuning_parameters_ct_equil_hpspevsv_node(app);
            create_tuning_parameters_ct_sd_node(app);
            create_tuning_parameters_ct_rocket_node(app);
        end
        
        function create_miscellaneous_node(app)
            % Create tuning_parameters_node
            app.miscellaneous_node = uitreenode(app.combustion_toolbox_node);
            app.miscellaneous_node.Text = 'Miscellaneous';
            app.miscellaneous_node.Tag = 'miscellaneous';
        end

        function create_tuning_parameters_ct_equil_tptv_node(app)
            % Create tuning_parameters_node subnode
            app.tuning_parameters_ct_equil_tptv_node = uitreenode(app.tuning_parameters_node);
            app.tuning_parameters_ct_equil_tptv_node.Text = 'CT-EQUIL: TP/TV';
            app.tuning_parameters_ct_equil_tptv_node.Tag = 'tuning_ct_equil_tptv';
        end
    
        function create_tuning_parameters_ct_equil_hpspevsv_node(app)
            % Create tuning_parameters_node subnode
            app.tuning_parameters_ct_equil_tptv_node = uitreenode(app.tuning_parameters_node);
            app.tuning_parameters_ct_equil_tptv_node.Text = 'CT-EQUIL: HP/SP/EV/SV';
            app.tuning_parameters_ct_equil_tptv_node.Tag = 'tuning_ct_equil_hpspevsv';
        end

        function create_tuning_parameters_ct_sd_node(app)
            % Create tuning_parameters_node subnode
            app.tuning_parameters_ct_sd_node = uitreenode(app.tuning_parameters_node);
            app.tuning_parameters_ct_sd_node.Text = 'CT-SD';
            app.tuning_parameters_ct_sd_node.Tag = 'tuning_ct_sd';
        end

        function create_tuning_parameters_ct_rocket_node(app)
            % Create tuning_parameters_node subnode
            app.tuning_parameters_ct_rocket_node = uitreenode(app.tuning_parameters_node);
            app.tuning_parameters_ct_rocket_node.Text = 'CT-ROCKET';
            app.tuning_parameters_ct_rocket_node.Tag = 'tuning_ct_rocket';
        end

        function fill_combustion_toolbox_node(app)
            % Fill combustion_toolbox_node

            % Create title1_label
           create_title(app, 'Combustion Toolbox Preferences')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Select an element in the tree to set its preferences.';
        end

        function fill_general_node(app)
            % Fill combustion_toolbox_node

            % Create title1_label
           create_title(app, 'Combustion Toolbox General Settings')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Set General settings of CT.';
        end

        function fill_tuning_parameters(app)
            % Fill combustion_toolbox_node

            % Create title1_label
           create_title(app, 'Combustion Toolbox Tuning Parameters')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the different CT modules.';
        end

        function fill_tuning_parameters_ct_equil_tptv_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Tuning Parameters: CT-EQUIL - TP/TV')

            % Tuning parameters: CT-EQUIL
            obj_num = obj_num + 1;
            obj_total = 5;
            % * Tuning parameters: CT-EQUIL - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-EQUIL - TP/TV', 'obj_num', obj_num);
            % ** Tuning parameters: CT-EQUIL - itMax_gibbs
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax_gibbs'};
            title = 'Maximum number of iterations of Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - itMax_ions
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax_ions'};
            title = 'Maximum number of iterations of charge balance (ions)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - tolN
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tolN'};
            title = 'Tolerance of the Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - tolE
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tolE'};
            title = 'Tolerance of the mass balance';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-EQUIL - tol_pi_e
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_pi_e'};
            title = 'Tolerance of the dimensionless Lagrangian multiplier - ions';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
        end
        
        function fill_tuning_parameters_ct_equil_hpspevsv_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Tuning Parameters: CT-EQUIL - HP/SP/EV/SV')

            % Tuning parameters: CT-EQUIL
            obj_num = obj_num + 1;
            obj_total = 6;
            % * Tuning parameters: CT-EQUIL - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-EQUIL - HP/SP/EV/SV', 'obj_num', obj_num);

            % ** Tuning parameters: CT-EQUIL - tol0
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol0'};
            title = 'Tolerance of the root finding algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - root_method
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_method'};
            title = 'Root finding method';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'newton', 'steff'}, 'subtype', 'function')

            % ** Tuning parameters: CT-EQUIL - itMax
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax'};
            title = 'Maximum number of iterations for root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-EQUIL - root_T0_l
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0_l'};
            title = 'First guess T[K] left branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%d K')

            % ** Tuning parameters: CT-EQUIL - root_T0_r
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0_r'};
            title = 'First guess T[K] right branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%d K')

            % ** Tuning parameters: CT-EQUIL - root_T0
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0'};
            title = 'Guess T[K] if it''s of previous range - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%d K')
            
        end

        function fill_tuning_parameters_ct_sd_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Tuning Parameters: CT-SD')

            % Tuning parameters: CT-SD
            obj_num = obj_num + 1;
            obj_total = 7;
            % * Tuning parameters: CT-SD - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-SD', 'obj_num', obj_num);
            % ** Tuning parameters: CT-SD - tol_shocks
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_shocks'};
            title = 'Tolerance of shocks/detonations kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-SD - it_shocks
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'it_shocks'};
            title = 'Maximum number of iterations for shocks and detonations kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-SD - Mach_thermo
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'Mach_thermo'};
            title = 'Preshock Mach number above which T2_guess will be computed from Eq.(XX)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.2g')
            
            % ** Tuning parameters: CT-SD - tol_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_oblique'};
            title = 'Tolerance oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'it_oblique'};
            title = 'Maximum number of iterations for oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - N_points_polar
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'N_points_polar'};
            title = 'Number of points to compute shock/detonation polar curves';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_guess_det
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'it_guess_det'};
            title = 'Maximum number of iterations to calculate detonation guess values';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
        end

        function fill_tuning_parameters_ct_rocket_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Tuning Parameters: CT-ROCKET')

            % Tuning parameters: CT-SD
            obj_num = obj_num + 1;
            obj_total = 2;
            % * Tuning parameters: CT-SD - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-ROCKET', 'obj_num', obj_num);
            % ** Tuning parameters: CT-ROCKET - tol_rocket
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_rocket'};
            title = 'Tolerance of rocket algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-ROCKET - it_rocket
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'it_rocket'};
            title = 'Maximum number of iterations for rocket algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
        end

        function fill_miscellaneous_node(app)
            % Fill combustion_toolbox_node

            % Create title1_label
           create_title(app, 'Combustion Toolbox Miscellaneous parameters')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Set miscellaneous parameters of CT.';
        end

        function create_ok_button(app)
            % Create ok_button
            app.ok_button = uibutton(app.preferences_UIFigure, 'push');
            app.ok_button.Position = [437 11 90 22];
            app.ok_button.Text = 'OK';
        end

        function create_cancel_button(app)
            % Create cancel_button
            app.cancel_button = uibutton(app.preferences_UIFigure, 'push');
            app.cancel_button.Position = [537 11 90 22];
            app.cancel_button.Text = 'Cancel';
            app.cancel_button.ButtonPushedFcn = createCallbackFcn(app, @cancel_button_pushed, true);
        end

        function create_apply_button(app)
            % Create apply_button
            app.apply_button = uibutton(app.preferences_UIFigure, 'push');
            app.apply_button.Position = [637 11 90 22];
            app.apply_button.Text = 'Apply';
        end

    end

end
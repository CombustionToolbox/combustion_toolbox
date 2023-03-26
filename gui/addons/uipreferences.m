classdef uipreferences < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        preferences_UIFigure                     matlab.ui.Figure
        cancel_button                            matlab.ui.control.Button
        ok_button                                matlab.ui.control.Button
        label_version                            matlab.ui.control.Label
        tree                                     matlab.ui.container.Tree
        combustion_toolbox_node                  matlab.ui.container.TreeNode
        general_node                             matlab.ui.container.TreeNode
        constants_node                           matlab.ui.container.TreeNode
        tuning_parameters_node                   matlab.ui.container.TreeNode
        tuning_parameters_flags_node             matlab.ui.container.TreeNode
        tuning_parameters_ct_equil_tptv_node     matlab.ui.container.TreeNode
        tuning_parameters_ct_equil_hpspevsv_node matlab.ui.container.TreeNode
        tuning_parameters_ct_sd_node             matlab.ui.container.TreeNode
        tuning_parameters_ct_rocket_node         matlab.ui.container.TreeNode
        miscellaneous_node                       matlab.ui.container.TreeNode
        miscellaneous_flags_node                 matlab.ui.container.TreeNode
        miscellaneous_plots_node                 matlab.ui.container.TreeNode
        miscellaneous_axes_node                  matlab.ui.container.TreeNode
        miscellaneous_export_node                matlab.ui.container.TreeNode
        ContextMenu                              matlab.ui.container.ContextMenu
        SnapshotMenu                             matlab.ui.container.Menu
    end
    
    properties (Access = public)
        C % Constants
        TN % Tuning properties
        Misc % Miscellaneous properties
        x0_panel_right = 206 % Initial positition of right panel in the x-axis [pixels]
        y0_panel_right = 428 % Initial positition of right panel in the y-axis [pixels]
        height_panel_0 = 38 % Default height of the right panel [pixels]
        delta_x = 9 % Left margin in the right panel
        delta_y = -12 % Top margin in the right panel
        width_amplification = 1.15
        height_amplification = 1
        width_box = 60 % Default box width;
        height_box = 22 % Default box height;
        height_text = 14 % Default box height;
        width_right = 521 % Default width right panel
        background_color = [0.9098 0.9098 0.8902] % Backgound color of the app
        caller_app % Handle to caller app
        dynamic_components % Struct with all the dynamic components
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
            % Copy properties to avoid overheads
            app.C = app.caller_app.C;
            app.TN = app.caller_app.TN;
            app.Misc = app.caller_app.Misc;
            % Update label version
            app.label_version.Text = sprintf('Version: %s', app.C.release);
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
                case 'constants'
                    fill_constants_node(app);
                case 'tuning'
                    fill_tuning_parameters(app);
                case 'tuning_flags'
                    fill_tuning_parameters_flags_node(app);
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
                case 'miscellaneous_flags'
                    fill_miscellaneous_flags_node(app);
                case 'miscellaneous_plots'
                    fill_miscellaneous_plots_node(app);
                case 'miscellaneous_axes'
                    fill_miscellaneous_axes_node(app);
                case 'miscellaneous_export'
                    fill_miscellaneous_export_node(app);
            end
        end
        
        % Save changes and close app
        function ok_button_pushed(app, event)
            % Delete app
            delete(app)
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
        function set_value_inner(app, event, varargin)
            Tag = get_tag(app, event);
            app.caller_app.(Tag{1}).(Tag{2}).(Tag{3}) = event.Value;
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
            % Set horizontal alignment
            app.dynamic_components.(label_name).HorizontalAlignment = 'left';
            % Set position
            app.dynamic_components.(label_name).Position = [app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 width_label, height_label];
            % Set title
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
            
            % Create object
            app.dynamic_components.(edit_name) = uieditfield(app.dynamic_components.(panel_name), type);

            if strcmpi(type, 'numeric')
                
                if strcmpi(subtype, 'inner')
                    % Get value
                    value = app.(tag{1}).(tag{2}).(tag{3});
                    % Assign value changed function
                    app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value_inner, true);
                else
                    % Get value
                    value = app.(tag{1}).(tag{2});
                    % Assign value changed function
                    app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value, true);
                end
                
                % Set limits value
                app.dynamic_components.(edit_name).Limits = limits;
                % Set display format
                app.dynamic_components.(edit_name).ValueDisplayFormat = display_format;
                % Set value
                app.dynamic_components.(edit_name).Value = value;
                
            else
                switch lower(subtype)
                    case 'function'
                        % Set value
                        app.dynamic_components.(edit_name).Value = func2str(app.(tag{1}).(tag{2}));
                        % Assign value changed function
                        app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value_function, true);
                    case 'inner'
                        % Set value
                        app.dynamic_components.(edit_name).Value = app.(tag{1}).(tag{2}).(tag{3});
                        % Assign value changed function
                        app.dynamic_components.(edit_name).ValueChangedFcn = createCallbackFcn(app, @set_value_inner, true);
                end
            end
            
            app.dynamic_components.(edit_name).Position = [width_label + 2 * app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 app.width_box, app.height_box];
            % Set tag
            app.dynamic_components.(edit_name).Tag = fullfile(tag{1}, tag{2});
        end
        
        % Create dropdown panel's children
        function create_panel_dropdown(app, panel_name, obj_num, tag, title, items, varargin)
            % Create label panel's children
            width_label = create_panel_label(app, panel_name, obj_num, title);

            % Definitions
            dropdown_name = [panel_name, 'dropdown', sprintf('%d', obj_num)];
            % Default
            subtype = 'none';
            % Unpack
            for i = 1:2:nargin-6
                switch lower(varargin{i})
                    case 'subtype'
                        subtype = varargin{i+1};
                end
            end
            
            % Create object
            app.dynamic_components.(dropdown_name) = uidropdown(app.dynamic_components.(panel_name));
            % Set items
            app.dynamic_components.(dropdown_name).Items = items;
            % Set position
            app.dynamic_components.(dropdown_name).Position = [width_label + 2 * app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 1.2 * app.width_box, app.height_box];

            switch lower(subtype)
                case 'function'
                    % Set value
                    app.dynamic_components.(dropdown_name).Value = func2str(app.(tag{1}).(tag{2}));
                    % Assign value changed function
                    app.dynamic_components.(dropdown_name).ValueChangedFcn = createCallbackFcn(app, @set_value_function, true);
                case 'inner'
                    % Set value
                    app.dynamic_components.(dropdown_name).Value = app.(tag{1}).(tag{2}).(tag{3});
                    % Assign value changed function
                    app.dynamic_components.(dropdown_name).ValueChangedFcn = createCallbackFcn(app, @set_value_inner, true);
                    % Set tag
                    app.dynamic_components.(dropdown_name).Tag = fullfile(tag{1}, tag{2}, tag{3});
                    return

                otherwise
                    % Set value
                    app.dynamic_components.(dropdown_name).Value = app.(tag{1}).(tag{2});
                    % Assign value changed function
                    app.dynamic_components.(dropdown_name).ValueChangedFcn = createCallbackFcn(app, @set_value, true);
            end
            % Set tag
            app.dynamic_components.(dropdown_name).Tag = fullfile(tag{1}, tag{2});
        end

        % Create checkbox panel's children
        function create_panel_checkbox(app, panel_name, obj_num, tag, title)
            % Create label panel's children
            width_label = create_panel_label(app, panel_name, obj_num, title);

            % Definitions
            checkbox_name = [panel_name, 'checkbox', sprintf('%d', obj_num)];
            % Create object
            app.dynamic_components.(checkbox_name) = uicheckbox(app.dynamic_components.(panel_name), "Text", '');
            % Set value
            app.dynamic_components.(checkbox_name).Value = app.(tag{1}).(tag{2});
            % Set position
            app.dynamic_components.(checkbox_name).Position = [width_label + 2 * app.delta_x,...
                 app.dynamic_components.(panel_name).Position(4) + 3 * obj_num * app.delta_y,...
                 1.2 * app.width_box, app.height_box];
            % Assign value changed function
            app.dynamic_components.(checkbox_name).ValueChangedFcn = createCallbackFcn(app, @set_value, true);
            % Set tag
            app.dynamic_components.(checkbox_name).Tag = fullfile(tag{1}, tag{2});
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

        function SnapshotMenuSelected(app, event)
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
            
            % Create label version
            create_label_version(app);
            % Create tree
            create_tree(app);
            % Create combustion_toolbox_node
            create_combustion_toolbox_node(app);
            % Create general_node
            create_general_node(app);
            % Create constants_node
            create_constants_node(app);
            % Create tuning_parameters_node
            create_tuning_parameters_node(app);
            % Create miscellaneous_node
            create_miscellaneous_node(app);
            % Create ok_button
            create_ok_button(app);
            % Create cancel_button
            create_cancel_button(app);

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.preferences_UIFigure);

            % Create SnapshotMenu
            app.SnapshotMenu = uimenu(app.ContextMenu);
            app.SnapshotMenu.MenuSelectedFcn = createCallbackFcn(app, @SnapshotMenuSelected, true);
            app.SnapshotMenu.Text = 'Snapshot';

            % Assign app.ContextMenu
            app.preferences_UIFigure.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.preferences_UIFigure.Visible = 'on';
        end
        
        function create_label_version(app)
            % Create VersionLabel
            app.label_version = uilabel(app.preferences_UIFigure);
            app.label_version.Position = [12 11 178 22];
            app.label_version.Text = 'Version:';
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
        
        function create_constants_node(app)
            % Create constants_node
            app.constants_node = uitreenode(app.combustion_toolbox_node);
            app.constants_node.Text = 'Constants';
            app.constants_node.Tag = 'constants';
        end

        function create_tuning_parameters_node(app)
            % Create tuning_parameters_node
            app.tuning_parameters_node = uitreenode(app.combustion_toolbox_node);
            app.tuning_parameters_node.Text = 'Tuning parameters';
            app.tuning_parameters_node.Tag = 'tuning';
            % Create subnodes
            create_tuning_parameters_flags_node(app);
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
             % Create subnodes
            create_miscellaneous_flags_node(app);
            create_miscellaneous_plots_node(app);
            create_miscellaneous_axes_node(app);
            create_miscellaneous_export_node(app);
        end

        function create_tuning_parameters_flags_node(app)
            % Create tuning_parameters_node subnode
            app.tuning_parameters_flags_node = uitreenode(app.tuning_parameters_node);
            app.tuning_parameters_flags_node.Text = 'Flags';
            app.tuning_parameters_flags_node.Tag = 'tuning_flags';
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

        function create_miscellaneous_flags_node(app)
            % Create miscellaneous_flags_node subnode
            app.miscellaneous_flags_node = uitreenode(app.miscellaneous_node);
            app.miscellaneous_flags_node.Text = 'Flags';
            app.miscellaneous_flags_node.Tag = 'miscellaneous_flags';
        end

        function create_miscellaneous_plots_node(app)
            % Create miscellaneous_plots_node subnode
            app.miscellaneous_plots_node = uitreenode(app.miscellaneous_node);
            app.miscellaneous_plots_node.Text = 'Plots';
            app.miscellaneous_plots_node.Tag = 'miscellaneous_plots';
        end

        function create_miscellaneous_axes_node(app)
            % Create miscellaneous_axes_node subnode
            app.miscellaneous_axes_node = uitreenode(app.miscellaneous_node);
            app.miscellaneous_axes_node.Text = 'Axes';
            app.miscellaneous_axes_node.Tag = 'miscellaneous_axes';
        end

        function create_miscellaneous_export_node(app)
            % Create miscellaneous_flags_node subnode
            app.miscellaneous_export_node = uitreenode(app.miscellaneous_node);
            app.miscellaneous_export_node.Text = 'Export';
            app.miscellaneous_export_node.Tag = 'miscellaneous_export';
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
            % Fill general_node

            % Create title1_label
            create_title(app, 'Combustion Toolbox General Settings')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Set General settings of CT.';
        end

        function fill_constants_node(app)
            % Fill constants_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Create title1_label
            create_title(app, 'Combustion Toolbox Constants parameters')
            
            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 4;
            
            % * Constants
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Constants', 'obj_num', obj_num);

            % ** Constants: R0
            obj_num_children = obj_num_children + 1;
            tag = {'C', 'R0'};
            title = 'Universal gas constant';
            app.width_box = 2 * app.width_box;
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%f J/(K mol)')
            app.width_box = app.width_box / 2;

            % ** Constants: gravity
            obj_num_children = obj_num_children + 1;
            tag = {'C', 'gravity'};
            title = 'Standard acceleration due to gravity';
            app.width_box = 1.8 * app.width_box;
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%f m/s')
            app.width_box = app.width_box / 1.8;

            % ** Constants: mintol_display
            obj_num_children = obj_num_children + 1;
            tag = {'C', 'mintol_display'};
            title = 'Minimum tolerance display species';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%1.2e')

            % ** Constants: composition_units
            obj_num_children = obj_num_children + 1;
            tag = {'C', 'composition_units'};
            title = 'Set units of the composition';
            app.width_box = 1.5 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'mol', 'molar fraction', 'mass fraction'})
            app.width_box = app.width_box / 1.5;
        end

        function fill_tuning_parameters(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'Combustion Toolbox Tuning Parameters')
            
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the different CT modules.';
            
            % Create text2_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text2_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text2_label.FontWeight = 'bold';
            app.dynamic_components.text2_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text2_label.Text = '  CT-EQUIL';
            
            % Create text3_label
            nlines = 6;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;

            label_text = [...
                '  Computes the composition at the equilibrium of multicomponent gas mixtures that undergo ', newline,...
                '  canonical thermochemical transformations from an initial state (reactants), defined by its', newline,...
                '  initial composition, temperature, and pressure, to a final state (products), defined by a set', newline,...
                '  of chemical species (gaseous—included ions—or pure condensed phase) and two', newline,...
                '  thermodynamic state functions, such as enthalpy and pressure, e.g., for isobaric', newline,...
                '  combustion processes.'
                ];

            app.dynamic_components.text3_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text3_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text3_label.Text = label_text;
            
            % Create text4_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text4_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text4_label.FontWeight = 'bold';
            app.dynamic_components.text4_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text4_label.Text = '  CT-SD';
            
            % Create text5_label
            nlines = 1;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;
            
            label_text = '  Solves steady-state shock and detonation waves in either normal or oblique incidence.';

            app.dynamic_components.text5_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text5_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text5_label.Text = label_text;
            
            % Create text6_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text6_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text6_label.FontWeight = 'bold';
            app.dynamic_components.text6_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text6_label.Text = '  CT-ROCKET';
            
            % Create text7_label
            nlines = 2;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;
            
            label_text = ['  Calculates the theoretical performance of rocket engines under highly idealized', newline,...
                          '  conditions.'];

            app.dynamic_components.text7_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text7_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text7_label.Text = label_text;
        end
        
        function fill_tuning_parameters_flags_node(app)
            % Fill tuning_parameters_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Tuning Parameters: FLAGS')

            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 3;

            % * Tuning parameters: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);

            % ** Tuning parameters: Flags - FLAG_FAST
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'FLAG_FAST'};
            title = 'Use the composition of the previous calculation as guess';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_TCHEM_FROZEN
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'FLAG_TCHEM_FROZEN'};
            title = 'Consider thermochemically frozen gas';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_EXTRAPOLATE
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'FLAG_EXTRAPOLATE'};
            title = 'Allow linear extrapolation of the polynomials fits';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
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
            obj_total = 8;
            % * Tuning parameters: CT-EQUIL - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-EQUIL - TP/TV', 'obj_num', obj_num);
            
            % ** Tuning parameters: CT-EQUIL - tolN
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tolN'};
            title = 'Tolerance of the composition of the mixture';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - tol_gibbs
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_gibbs'};
            title = 'Tolerance of the Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-EQUIL - itMax_gibbs
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax_gibbs'};
            title = 'Maximum number of iterations of Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-EQUIL - tolN_guess
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tolN_guess'};
            title = 'Tolerance of the composition of the mixture (guess)';
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

            % ** Tuning parameters: CT-EQUIL - itMax_ions
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax_ions'};
            title = 'Maximum number of iterations of charge balance (ions)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - T_ions
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'T_ions'};
            title = 'Minimum temperature [K] to consider ionized species';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')
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
            
            % ** Tuning parameters: CT-EQUIL - itMax
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'itMax'};
            title = 'Maximum number of iterations for root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - root_method
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_method'};
            title = 'Root finding method';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'newton', 'steff', 'nsteff'}, 'subtype', 'function')
            
            % ** Tuning parameters: CT-EQUIL - root_T0_l
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0_l'};
            title = 'First temperature guess left branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')

            % ** Tuning parameters: CT-EQUIL - root_T0_r
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0_r'};
            title = 'First temperature guess right branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')

            % ** Tuning parameters: CT-EQUIL - root_T0
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'root_T0'};
            title = 'Temperature guess if it''s outside previous range - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K') 
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
            title = 'Pre-shock Mach number above which T2_guess will be computed from Eq.(XX)';
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
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - tol_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'tol_limitRR'};
            title = 'Tolerance limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'TN', 'it_limitRR'};
            title = 'Maximum number of iterations for limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

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
        
        function fill_miscellaneous_flags_node(app)
            % Fill miscellaneous_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Miscellaneous parameters: FLAGS')

            % Miscellaneous: Flags
            obj_num = obj_num + 1;
            obj_total = 1;
            
            % * Miscellaneous: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);

            % ** Miscellaneous: Flags - FLAG_FAST
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'FLAG_RESULTS'};
            title = 'Show results in the command window (desktop or standalone versions)';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
        end
        
        function fill_miscellaneous_plots_node(app)
            % Fill miscellaneous_plots_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Miscellaneous parameters: Plots')

            % Miscellaneous: Plots
            obj_num = obj_num + 1;
            obj_total = 8;
            
            % * Miscellaneous: Plots
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Plots', 'obj_num', obj_num);
            
            % ** Miscellaneous: Plots - linestyle
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'linestyle'};
            title = 'Line style';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'-', '--', ':', '-.'}, 'subtype', 'inner')

            % ** Miscellaneous: Plots - symbolstyle
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'symbolstyle'};
            title = 'Symbol style';
            list_symbols = {'o', '+', '*', '.', 'x', '_', '|', 'square', 'diamond', '^', 'v', '>', '<', 'pentragram', 'hexagram'};
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_symbols, 'subtype', 'inner')

            % ** Miscellaneous: Plots - linewidth
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'linewidth'};
            title = 'Linewidth';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0.1, Inf], 'display_format', '%g', 'subtype', 'inner')

            % ** Miscellaneous: Plots - fontsize
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'fontsize'};
            title = 'Fontsize';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [8, Inf], 'display_format', '%d', 'subtype', 'inner')

            % ** Miscellaneous: Plots - colorpalette
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'colorpalette'};
            title = 'Color palette';
            list_colorpalette = {'Seaborn', 'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Normal', 'Rainbow', 'Paired_custom', 'Custom'};
            app.width_box = 1.2 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_colorpalette, 'subtype', 'inner')
            app.width_box = app.width_box / 1.2;

            % ** Miscellaneous: Plots - colorpaletteLenght
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'colorpaletteLenght'};
            title = 'Maximum number of colors in the palette';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d', 'subtype', 'inner')

            % ** Miscellaneous: Plots - label_type
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'label_type'};
            title = 'Type of label';
            app.width_box = 1.2 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'short', 'medium', 'long'}, 'subtype', 'inner')
            app.width_box = app.width_box / 1.2;

            % ** Miscellaneous: Plots - legend_location
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'legend_location'};
            title = 'Legend location';
            list_legendlocation = {'best', 'bestoutside', 'northeastoutside', 'northwestoutside', 'southeastoutside', 'southwestoutside', 'northoutside', 'southoutside', 'eastoutside', 'westoutside', 'northeast', 'northwest', 'southeast', 'southwest', 'north', 'south', 'east', 'west'};
            app.width_box = 1.8 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_legendlocation, 'subtype', 'inner')
            app.width_box = app.width_box / 1.8;
        end

        function fill_miscellaneous_axes_node(app)
            % Fill miscellaneous_axes_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Miscellaneous parameters: Axes')

            % Miscellaneous: Axes
            obj_num = obj_num + 1;
            obj_total = 9;
            
            % * Miscellaneous: Axes
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Axes', 'obj_num', obj_num);

            % ** Miscellaneous: Axes - box
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'box'};
            title = 'Display axes outline';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - grid
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'grid'};
            title = 'Display or hide axes grid lines';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - hold
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'hold'};
            title = 'Retain current plot when adding new plots';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - axis_x
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'axis_x'};
            title = 'Set x-axis limits';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'auto', 'tight', 'padded', 'fill', 'equal', 'image', 'square', 'vis3d'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - axis_y
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'axis_y'};
            title = 'Set y-axis limits';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'auto', 'tight', 'padded', 'fill', 'equal', 'image', 'square', 'vis3d'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - xscale
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'xscale'};
            title = 'Set x-axis scale';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'linear', 'log'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - yscale
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'yscale'};
            title = 'Set y-axis scale';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'linear', 'log'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - xdir
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'xdir'};
            title = 'Set x-axis direction';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'normal', 'reverse'}, 'subtype', 'inner')

            % ** Miscellaneous: Axes - ydir
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'config', 'ydir'};
            title = 'Set y-axis direction';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'normal', 'reverse'}, 'subtype', 'inner')
        end

        function fill_miscellaneous_export_node(app)
            % Fill miscellaneous_export_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox Miscellaneous parameters: Export')

            % Miscellaneous: Export
            obj_num = obj_num + 1;
            obj_total = 2;
            
            % * Miscellaneous: Export
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Export', 'obj_num', obj_num);

            % ** Miscellaneous: Export - format
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'export_results', 'format'};
            title = 'Default file format to export results';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'.xls', '.mat'}, 'subtype', 'inner')

            % ** Miscellaneous: Export - filename
            obj_num_children = obj_num_children + 1;
            tag = {'Misc', 'export_results', 'filename'};
            title = 'Default file name to export results';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'type', 'text', 'subtype', 'inner')
        end

        function create_ok_button(app)
            % Create ok_button
            app.ok_button = uibutton(app.preferences_UIFigure, 'push');
            app.ok_button.Position = [537 11 90 22];
            app.ok_button.Text = 'OK';
            app.ok_button.ButtonPushedFcn = createCallbackFcn(app, @ok_button_pushed, true);
        end

        function create_cancel_button(app)
            % Create cancel_button
            app.cancel_button = uibutton(app.preferences_UIFigure, 'push');
            app.cancel_button.Position = [637 11 90 22];
            app.cancel_button.Text = 'Cancel';
            app.cancel_button.ButtonPushedFcn = createCallbackFcn(app, @cancel_button_pushed, true);
        end

    end

end
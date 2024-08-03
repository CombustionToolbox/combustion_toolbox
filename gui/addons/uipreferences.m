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
        solvers_node                             matlab.ui.container.TreeNode
        equilibriumSolver_node                   matlab.ui.container.TreeNode
        equilibriumSolver_flags_node             matlab.ui.container.TreeNode
        equilibriumSolver_module1_node           matlab.ui.container.TreeNode
        equilibriumSolver_module2_node           matlab.ui.container.TreeNode
        shockSolver_node                         matlab.ui.container.TreeNode
        shockSolver_flags_node                   matlab.ui.container.TreeNode
        shockSolver_module1_node                 matlab.ui.container.TreeNode
        detonationSolver_node                    matlab.ui.container.TreeNode
        detonationSolver_flags_node              matlab.ui.container.TreeNode
        detonationSolver_module1_node            matlab.ui.container.TreeNode
        rocketSolver_node                        matlab.ui.container.TreeNode
        rocketSolver_module1_node                matlab.ui.container.TreeNode
        rocketSolver_flags_node                  matlab.ui.container.TreeNode
        plotConfig_node                          matlab.ui.container.TreeNode
        plotConfig_plots_node                 matlab.ui.container.TreeNode
        plotConfig_axes_node                  matlab.ui.container.TreeNode
        miscellaneous_export_node                matlab.ui.container.TreeNode
        ContextMenu                              matlab.ui.container.ContextMenu
        SnapshotMenu                             matlab.ui.container.Menu
    end
    
    properties (Access = public)
        constants            % Constants object
        equilibriumSolver    % EquilibriumSolver object
        shockSolver          % ShockSolver object
        detonationSolver     % DetonationSolver object
        rocketSolver         % RocketSolver object
        plotConfig           % PlotConfig object
        export               % Export object
        x0_panel_right = 206 % Initial positition of right panel in the x-axis [pixels]
        y0_panel_right = 428 % Initial positition of right panel in the y-axis [pixels]
        height_panel_0 = 38 % Default height of the right panel [pixels]
        delta_x = 9 % Left margin in the right panel
        delta_y = -12 % Top margin in the right panel
        width_amplification = 1.2
        height_amplification = 1
        width_box = 60 % Default box width;
        height_box = 22 % Default box height;
        height_text = 16 % Default box height;
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
                app.init;
            end

            % Copy properties to avoid overheads
            app.constants = app.caller_app.constants;
            app.equilibriumSolver = app.caller_app.equilibriumSolver;
            app.shockSolver = app.caller_app.shockSolver;
            app.detonationSolver = app.caller_app.detonationSolver;
            app.rocketSolver = app.caller_app.rocketSolver;
            app.plotConfig = app.caller_app.plotConfig;
            app.export = app.caller_app.export;

            % Update label version
            app.label_version.Text = sprintf('Version: %s', app.constants.release);

            % Create right panel based on first node uitree
            fill_combustion_toolbox_node(app);

            % Expand uitree
            expand(app.tree);
            
            % Update width_amplification (only macOS)
            % if ismac
            %     app.width_amplification = 0.91;
            % end
            
        end

        function init(app)
            % Import packages
            import combustiontoolbox.common.Constants
            import combustiontoolbox.equilibrium.*
            import combustiontoolbox.shockdetonation.*
            import combustiontoolbox.rocketSolver.*
            import combustiontoolbox.utils.display.PlotConfig
            import combustiontoolbox.utils.Export

            % Assign default values
            app.caller_app.constants = Constants();
            app.caller_app.plotConfig = PlotConfig();
            app.caller_app.export = Export();

            % Initialize solvers
            app.caller_app.equilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('plotConfig', app.caller_app.plotConfig);
            app.caller_app.shockSolver = combustiontoolbox.shockdetonation.ShockSolver('plotConfig', app.caller_app.plotConfig, 'equilibriumSolver', app.caller_app.equilibriumSolver);
            app.caller_app.detonationSolver = combustiontoolbox.shockdetonation.DetonationSolver('plotConfig', app.caller_app.plotConfig, 'equilibriumSolver', app.caller_app.equilibriumSolver);
            app.caller_app.rocketSolver = combustiontoolbox.rocket.RocketSolver('plotConfig', app.caller_app.plotConfig, 'equilibriumSolver', app.caller_app.equilibriumSolver);
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
                case 'solvers'
                    fill_solvers_node(app);
                case 'equilibriumsolver'
                    fill_equilibriumSolver_node(app);
                case 'equilibriumsolver_flags'
                    fill_equilibriumSolver_flags_node(app);
                case 'equilibriumsolver_module1'
                    fill_equilibriumSolver_module1_node(app);
                case 'equilibriumsolver_module2'
                    fill_equilibriumSolver_module2_node(app);
                case 'shocksolver'
                    fill_shockSolver_node(app);
                case 'shocksolver_flags'
                    fill_shockSolver_flags_node(app);
                case 'shocksolver_module1'
                    fill_shockSolver_module1_node(app);
                case 'detonationsolver'
                    fill_detonationSolver_node(app);
                case 'detonationsolver_flags'
                    fill_detonationSolver_flags_node(app);
                case 'detonationsolver_module1'
                    fill_detonationSolver_module1_node(app);
                case 'rocketsolver'
                    fill_rocketSolver_node(app);
                case 'rocketsolver_flags'
                    fill_rocketSolver_flags_node(app);
                case 'rocketsolver_module1'
                    fill_rocketSolver_module1_node(app);
                case 'plotconfig'
                    fill_plotConfig_node(app);
                case 'plotconfig_plots'
                    fill_plotConfig_plots_node(app);
                case 'plotconfig_axes'
                    fill_plotConfig_axes_node(app);
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
            init(app);

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
            % Create solvers_node
            create_solvers_node(app);
            % Create equilibriumSolver_node
            create_equilibriumSolver_node(app);
            % Create shockSolver_node
            create_shockSolver_node(app);
            % Create detonationSolver_node
            create_detonationSolver_node(app);
            % Create rocketSolver_node
            create_rocketSolver_node(app);
            % Create plotConfig_node
            create_plotConfig_node(app);
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

        function create_solvers_node(app)
            % Create solvers_node
            app.solvers_node = uitreenode(app.combustion_toolbox_node);
            app.solvers_node.Text = 'Solvers';
            app.solvers_node.Tag = 'solvers';
        end

        function create_equilibriumSolver_node(app)
            % Create equilibriumSolver_node
            app.equilibriumSolver_node = uitreenode(app.solvers_node);
            app.equilibriumSolver_node.Text = 'EquilibriumSolver';
            app.equilibriumSolver_node.Tag = 'equilibriumSolver';
            % Create subnodes
            create_equilibriumSolver_flags_node(app);
            create_equilibriumSolver_module1_node(app);
            create_equilibriumSolver_module2_node(app);
        end

        function create_equilibriumSolver_flags_node(app)
            % Create tuning_parameters_node subnode
            app.equilibriumSolver_flags_node = uitreenode(app.equilibriumSolver_node);
            app.equilibriumSolver_flags_node.Text = 'Flags';
            app.equilibriumSolver_flags_node.Tag = 'equilibriumSolver_flags';
        end

        function create_equilibriumSolver_module1_node(app)
            % Create tuning_parameters_node subnode
            app.equilibriumSolver_module1_node = uitreenode(app.equilibriumSolver_node);
            app.equilibriumSolver_module1_node.Text = 'TP/TV problems';
            app.equilibriumSolver_module1_node.Tag = 'equilibriumSolver_module1';
        end
    
        function create_equilibriumSolver_module2_node(app)
            % Create tuning_parameters_node subnode
            app.equilibriumSolver_module1_node = uitreenode(app.equilibriumSolver_node);
            app.equilibriumSolver_module1_node.Text = 'HP/SP/EV/SV problems';
            app.equilibriumSolver_module1_node.Tag = 'equilibriumSolver_module2';
        end
        
        function create_shockSolver_node(app)
            % Create equilibriumSolver_node
            app.shockSolver_node = uitreenode(app.solvers_node);
            app.shockSolver_node.Text = 'ShockSolver';
            app.shockSolver_node.Tag = 'shockSolver';
            % Create subnodes
            create_shockSolver_flags_node(app);
            create_shockSolver_module1_node(app);
        end

        function create_shockSolver_flags_node(app)
            % Create tuning_parameters_node subnode
            app.shockSolver_flags_node = uitreenode(app.shockSolver_node);
            app.shockSolver_flags_node.Text = 'Flags';
            app.shockSolver_flags_node.Tag = 'shockSolver_flags';
        end

        function create_shockSolver_module1_node(app)
            % Create tuning_parameters_node subnode
            app.shockSolver_module1_node = uitreenode(app.shockSolver_node);
            app.shockSolver_module1_node.Text = 'Shock problems';
            app.shockSolver_module1_node.Tag = 'shockSolver_module1';
        end

        function create_detonationSolver_node(app)
            % Create equilibriumSolver_node
            app.detonationSolver_node = uitreenode(app.solvers_node);
            app.detonationSolver_node.Text = 'DetonationSolver';
            app.detonationSolver_node.Tag = 'detonationSolver';
            % Create subnodes
            create_detonationSolver_flags_node(app);
            create_detonationSolver_module1_node(app);
        end

        function create_detonationSolver_flags_node(app)
            % Create tuning_parameters_node subnode
            app.detonationSolver_flags_node = uitreenode(app.detonationSolver_node);
            app.detonationSolver_flags_node.Text = 'Flags';
            app.detonationSolver_flags_node.Tag = 'detonationSolver_flags';
        end

        function create_detonationSolver_module1_node(app)
            % Create tuning_parameters_node subnode
            app.detonationSolver_module1_node = uitreenode(app.detonationSolver_node);
            app.detonationSolver_module1_node.Text = 'Detonation problems';
            app.detonationSolver_module1_node.Tag = 'detonationSolver_module1';
        end

        function create_rocketSolver_node(app)
            % Create equilibriumSolver_node
            app.rocketSolver_node = uitreenode(app.solvers_node);
            app.rocketSolver_node.Text = 'RocketSolver';
            app.rocketSolver_node.Tag = 'rocketSolver';
            % Create subnodes
            create_rocketSolver_flags_node(app);
            create_rocketSolver_module1_node(app);
        end

        function create_rocketSolver_flags_node(app)
            % Create tuning_parameters_node subnode
            app.rocketSolver_flags_node = uitreenode(app.rocketSolver_node);
            app.rocketSolver_flags_node.Text = 'Flags';
            app.rocketSolver_flags_node.Tag = 'rocketSolver_flags';
        end

        function create_rocketSolver_module1_node(app)
            % Create tuning_parameters_node subnode
            app.rocketSolver_module1_node = uitreenode(app.rocketSolver_node);
            app.rocketSolver_module1_node.Text = 'Rocket problems';
            app.rocketSolver_module1_node.Tag = 'rocketSolver_module1';
        end

        function create_plotConfig_node(app)
            % Create plotConfig_node
            app.plotConfig_node = uitreenode(app.combustion_toolbox_node);
            app.plotConfig_node.Text = 'Plot configuration';
            app.plotConfig_node.Tag = 'plotConfig';
             % Create subnodes
            create_plotConfig_plots_node(app);
            create_plotConfig_axes_node(app);
        end

        function create_plotConfig_plots_node(app)
            % Create plotConfig_plots_node subnode
            app.plotConfig_plots_node = uitreenode(app.plotConfig_node);
            app.plotConfig_plots_node.Text = 'Plots';
            app.plotConfig_plots_node.Tag = 'plotConfig_plots';
        end

        function create_plotConfig_axes_node(app)
            % Create plotConfig_axes_node subnode
            app.plotConfig_axes_node = uitreenode(app.plotConfig_node);
            app.plotConfig_axes_node.Text = 'Axes';
            app.plotConfig_axes_node.Tag = 'plotConfig_axes';
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
            app.dynamic_components.text1_label.Text = '  Set general settings of CT.';
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
            obj_total = 3;
            
            % * Constants
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Constants', 'obj_num', obj_num);

            % ** Constants: R0
            obj_num_children = obj_num_children + 1;
            tag = {'constants', 'R0'};
            title = 'Universal gas constant';
            app.width_box = 2 * app.width_box;
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%f J/(K mol)')
            app.width_box = app.width_box / 2;

            % ** Constants: gravity
            obj_num_children = obj_num_children + 1;
            tag = {'constants', 'G'};
            title = 'Standard acceleration due to gravity';
            app.width_box = 1.8 * app.width_box;
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%f m/s2')
            app.width_box = app.width_box / 1.8;

            % ** Constants: gravity
            obj_num_children = obj_num_children + 1;
            tag = {'constants', 'NA'};
            title = 'Avogadro constant';
            app.width_box = 2.7 * app.width_box;
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.4e molecule/mol')
            app.width_box = app.width_box / 2.7;

            % % ** Constants: mintol_display
            % obj_num_children = obj_num_children + 1;
            % tag = {'constants', 'mintolDisplay'};
            % title = 'Minimum tolerance display species';
            % create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%1.2e')
            % 
            % % ** Constants: composition_units
            % obj_num_children = obj_num_children + 1;
            % tag = {'constants', 'compositionUnits'};
            % title = 'Set units of the composition';
            % app.width_box = 1.5 * app.width_box;
            % create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'mol', 'molar fraction', 'mass fraction'})
            % app.width_box = app.width_box / 1.5;
        end
        
        function fill_solvers_node(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'Solvers components settings')
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the different CT modules.';
        end

        function fill_equilibriumSolver_node(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'EquilibriumSolver class')
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the EquilibriumSolver class.';
            
            % Create text2_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text2_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text2_label.FontWeight = 'bold';
            app.dynamic_components.text2_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text2_label.Text = '  CT-EQUIL';
            
            % Create text3_label
            nlines = 8;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;

            label_text = [
                '  The EquilibriumSolver class, referred to as the CT-EQUIL module, serves as the foundational', newline,...
                '  computational engine for all solvers within the system. It calculates the equilibrium', newline,...
                '  composition of multicomponent gas mixtures undergoing canonical thermochemical', newline,...
                '  transformations. These transformations occur from an initial state (reactants), defined by', newline,...
                '  its initial composition, temperature, and pressure, to a final state (products), characterized', newline,...
                '  by a set of chemical species (gaseous, including ions, or pure condensed phase) and two', newline,...
                '  thermodynamic state functions, such as enthalpy and pressure. This is applicable to processes', newline,...
                '  like isobaric combustion.'
            ];

            app.dynamic_components.text3_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text3_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text3_label.Text = label_text;
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
        
        function fill_equilibriumSolver_flags_node(app)
            % Fill tuning_parameters_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'EquilibriumSolver class: FLAGS')

            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 4;

            % * Tuning parameters: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);
            
            % ** Problem Description: Flags - FLAG_TCHEM_FROZEN
            % obj_num_children = obj_num_children + 1;
            % tag = {'equilibriumSolver', 'FLAG_TCHEM_FROZEN'};
            % title = 'Consider thermochemically frozen gas';
            % create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_FAST
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'FLAG_FAST'};
            title = 'Use the composition of the previous calculation as guess';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_EXTRAPOLATE
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'FLAG_EXTRAPOLATE'};
            title = 'Allow linear extrapolation of the polynomials fits';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_RESULTS
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'FLAG_RESULTS'};
            title = 'Print results in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_TIME
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'FLAG_TIME'};
            title = 'Print elapsed time in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
        end

        function fill_equilibriumSolver_module1_node(app)
            % Fill tuning_parameters_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'EquilibriumSolver class: TP/TV')

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
            tag = {'equilibriumSolver', 'tolMoles'};
            title = 'Tolerance of the composition of the mixture';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - tol_gibbs
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'tolGibbs'};
            title = 'Tolerance of the Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-EQUIL - itMax_gibbs
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'itMaxGibbs'};
            title = 'Maximum number of iterations of Gibbs/Helmholtz minimization method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-EQUIL - tolN_guess
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'tolMolesGuess'};
            title = 'Tolerance of the composition of the mixture (guess)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - tolE
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'tolE'};
            title = 'Tolerance of the mass balance';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - tol_pi_e
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'tolMultiplierIons'};
            title = 'Tolerance of the dimensionless Lagrangian multiplier - ions';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-EQUIL - itMax_ions
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'itMaxIons'};
            title = 'Maximum number of iterations of charge balance (ions)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - T_ions
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'temperatureIons'};
            title = 'Minimum temperature [K] to consider ionized species';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')
        end
        
        function fill_equilibriumSolver_module2_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'EquilibriumSolver class: HP/SP/EV/SV')

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
            tag = {'equilibriumSolver', 'tol0'};
            title = 'Tolerance of the root finding algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-EQUIL - itMax
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'itMax'};
            title = 'Maximum number of iterations for root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-EQUIL - root_method
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'rootMethod'};
            title = 'Root finding method';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'newton', 'steff', 'nsteff'}, 'subtype', 'function')
            
            % ** Tuning parameters: CT-EQUIL - root_T0_l
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'root_T0_l'};
            title = 'First temperature guess left branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')

            % ** Tuning parameters: CT-EQUIL - root_T0_r
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'root_T0_r'};
            title = 'First temperature guess right branch - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K')

            % ** Tuning parameters: CT-EQUIL - root_T0
            obj_num_children = obj_num_children + 1;
            tag = {'equilibriumSolver', 'root_T0'};
            title = 'Temperature guess if it''s outside previous range - root finding method';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0, Inf], 'display_format', '%g K') 
        end

        function fill_shockSolver_node(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'ShockSolver class')
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the ShockSolver class.';
            
            % Create text2_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text2_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text2_label.FontWeight = 'bold';
            app.dynamic_components.text2_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text2_label.Text = '  CT-SD';
            
            % Create text3_label
            nlines = 2;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;

            label_text = [
                '  The ShockSolver class, from the CT-SD module, solves steady-state shock', newline,...
                '  either normal or oblique incidence.',...
            ];

            app.dynamic_components.text3_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text3_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text3_label.Text = label_text;
        end

        function fill_shockSolver_flags_node(app)
            % Fill tuning_parameters_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'ShockSolver class: FLAGS')

            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 2;

            % * Tuning parameters: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);
            
            % ** Problem Description: Flags - FLAG_TCHEM_FROZEN
            % obj_num_children = obj_num_children + 1;
            % tag = {'equilibriumSolver', 'FLAG_TCHEM_FROZEN'};
            % title = 'Consider thermochemically frozen gas';
            % create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_RESULTS
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'FLAG_RESULTS'};
            title = 'Print results in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_TIME
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'FLAG_TIME'};
            title = 'Print elapsed time in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
        end

        function fill_shockSolver_module1_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'ShockSolver class')

            % Tuning parameters: CT-SD
            obj_num = obj_num + 1;
            obj_total = 6;

            % * Tuning parameters: CT-SD - panel
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Module: CT-SD', 'obj_num', obj_num);

            % ** Tuning parameters: CT-SD - tol_shocks
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'tol0'};
            title = 'Tolerance of shock kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-SD - it_shocks
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'itMax'};
            title = 'Maximum number of iterations for shocks and detonations kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-SD - Mach_thermo
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'machThermo'};
            title = 'Pre-shock Mach number above which T2_guess will be computed from Eq.(XX)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.2g')
            
            % ** Tuning parameters: CT-SD - tol_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'tolOblique'};
            title = 'Tolerance oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'itOblique'};
            title = 'Maximum number of iterations for oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - N_points_polar
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'numPointsPolar'};
            title = 'Number of points to compute shock/detonation polar curves';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - tol_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'tolLimitRR'};
            title = 'Tolerance limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'shockSolver', 'itLimitRR'};
            title = 'Maximum number of iterations for limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
        end

        function fill_detonationSolver_node(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'DetonationSolver class')
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the DetonationSolver class.';
            
            % Create text2_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text2_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text2_label.FontWeight = 'bold';
            app.dynamic_components.text2_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text2_label.Text = '  CT-SD';
            
            % Create text3_label
            nlines = 3;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;

            label_text = [
                '  The DetonationSolver class, from the CT-SD module, solves steady-state detonations', newline,...
                '  either normal or oblique incidence. Overdriven and underdriven detonations are also', newline,...
                '  considered.'
            ];

            app.dynamic_components.text3_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text3_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text3_label.Text = label_text;
        end

        function fill_detonationSolver_flags_node(app)
            % Fill tuning_parameters_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'DetonationSolver class: FLAGS')

            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 2;

            % * Tuning parameters: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);

            % ** Tuning parameters: Flags - FLAG_RESULTS
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'FLAG_RESULTS'};
            title = 'Print results in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_TIME
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'FLAG_TIME'};
            title = 'Print elapsed time in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
        end

        function fill_detonationSolver_module1_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'DetonationSolver class')

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
            tag = {'detonationSolver', 'tol0'};
            title = 'Tolerance of shock kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-SD - it_shocks
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'itMax'};
            title = 'Maximum number of iterations for shocks and detonations kernel';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-SD - Mach_thermo
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'machThermo'};
            title = 'Pre-shock Mach number above which T2_guess will be computed from Eq.(XX)';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.2g')
            
            % ** Tuning parameters: CT-SD - tol_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'tolOblique'};
            title = 'Tolerance oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_oblique
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'itOblique'};
            title = 'Maximum number of iterations for oblique shocks algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - N_points_polar
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'numPointsPolar'};
            title = 'Number of points to compute shock/detonation polar curves';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
            
            % ** Tuning parameters: CT-SD - tol_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'tolLimitRR'};
            title = 'Tolerance limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')
            
            % ** Tuning parameters: CT-SD - it_limitRR
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'itLimitRR'};
            title = 'Maximum number of iterations for limit of regular reflections algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Tuning parameters: CT-SD - itGuess
            obj_num_children = obj_num_children + 1;
            tag = {'detonationSolver', 'itGuess'};
            title = 'Maximum number of iterations to calculate detonation guess values';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
        end

        function fill_rocketSolver_node(app)
            % Fill tuning_parameters_node

            % Create title1_label
            create_title(app, 'RocketSolver class')
            
            % Create text1_label
            y0 = 417;
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text1_label.Text = '  Set tuning parameters for the RocketSolver class.';
            
            % Create text2_label
            nlines = 1.2;
            y0 = y0 - app.height_box * nlines;
            app.dynamic_components.text2_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text2_label.FontWeight = 'bold';
            app.dynamic_components.text2_label.Position = [206, y0, app.width_right, 22];
            app.dynamic_components.text2_label.Text = '  CT-ROCKET';
            
            % Create text3_label
            nlines = 2;
            y0 = y0 - app.height_text * nlines;
            h = app.height_text * nlines + 1;

            label_text = [
                '  The RocketSolver class, from the CT-Rocket module, calculates the theoretical', newline,...
                '  performance of rocket engines under highly idealized conditions.'
            ];

            app.dynamic_components.text3_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text3_label.Position = [206, y0, app.width_right, h];
            app.dynamic_components.text3_label.Text = label_text;
        end

        function fill_rocketSolver_flags_node(app)
            % Fill tuning_parameters_flags_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'RocketSolver class: FLAGS')

            % Tuning parameters: Flags
            obj_num = obj_num + 1;
            obj_total = 2;

            % * Tuning parameters: Flags
            panel_name = create_panel(app, [app.x0_panel_right,...
                app.y0_panel_right - obj_total * app.height_panel_0,...
                app.width_right,  obj_total * app.height_panel_0],...
                ' Flags', 'obj_num', obj_num);
            
            % ** Problem Description: Flags - FLAG_TCHEM_FROZEN
            % obj_num_children = obj_num_children + 1;
            % tag = {'equilibriumSolver', 'FLAG_TCHEM_FROZEN'};
            % title = 'Consider thermochemically frozen gas';
            % create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_RESULTS
            obj_num_children = obj_num_children + 1;
            tag = {'rocketSolver', 'FLAG_RESULTS'};
            title = 'Print results in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)

            % ** Tuning parameters: Flags - FLAG_TIME
            obj_num_children = obj_num_children + 1;
            tag = {'rocketSolver', 'FLAG_TIME'};
            title = 'Print elapsed time in the command window';
            create_panel_checkbox(app, panel_name, obj_num_children, tag, title)
        end
        
        function fill_rocketSolver_module1_node(app)
            % Fill tuning_parameters_node
            
            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'RocketSolver class')

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
            tag = {'rocketSolver', 'tol0'};
            title = 'Tolerance of rocket algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1e-100, Inf], 'display_format', '%.1e')

            % ** Tuning parameters: CT-ROCKET - it_rocket
            obj_num_children = obj_num_children + 1;
            tag = {'rocketSolver', 'itMax'};
            title = 'Maximum number of iterations for rocket algorithm';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')
        end

        function fill_plotConfig_node(app)
            % Fill combustion_toolbox_node

            % Create title1_label
           create_title(app, 'Combustion Toolbox plot configuration')

            % Create text1_label
            app.dynamic_components.text1_label = uilabel(app.preferences_UIFigure);
            app.dynamic_components.text1_label.Position = [206 417 app.width_right 22];
            app.dynamic_components.text1_label.Text = '  Set the configuration of the plots of CT.';
        end
        
        function fill_plotConfig_plots_node(app)
            % Fill plotConfig_plots_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox plot configuration: Plots')

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
            tag = {'plotConfig', 'linestyle'};
            title = 'Line style';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'-', '--', ':', '-.'})

            % ** Miscellaneous: Plots - symbolstyle
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'symbolstyle'};
            title = 'Symbol style';
            list_symbols = {'o', '+', '*', '.', 'x', '_', '|', 'square', 'diamond', '^', 'v', '>', '<', 'pentragram', 'hexagram'};
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_symbols)

            % ** Miscellaneous: Plots - linewidth
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'linewidth'};
            title = 'Linewidth';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [0.1, Inf], 'display_format', '%g')

            % ** Miscellaneous: Plots - fontsize
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'fontsize'};
            title = 'Fontsize';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [8, Inf], 'display_format', '%d')

            % ** Miscellaneous: Plots - colorpalette
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'colorpalette'};
            title = 'Color palette';
            list_colorpalette = {'Seaborn', 'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3', 'Blues', 'BuGn', 'BuPu', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'Normal', 'Rainbow', 'Paired_custom', 'Custom'};
            app.width_box = 1.2 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_colorpalette)
            app.width_box = app.width_box / 1.2;

            % ** Miscellaneous: Plots - colorpaletteLenght
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'colorpaletteLenght'};
            title = 'Maximum number of colors in the palette';
            create_panel_edit_field(app, panel_name, obj_num_children, tag, title, 'limits', [1, Inf], 'display_format', '%d')

            % ** Miscellaneous: Plots - label_type
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'label_type'};
            title = 'Type of label';
            app.width_box = 1.2 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'short', 'medium', 'long'})
            app.width_box = app.width_box / 1.2;

            % ** Miscellaneous: Plots - legend_location
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'legend_location'};
            title = 'Legend location';
            list_legendlocation = {'best', 'bestoutside', 'northeastoutside', 'northwestoutside', 'southeastoutside', 'southwestoutside', 'northoutside', 'southoutside', 'eastoutside', 'westoutside', 'northeast', 'northwest', 'southeast', 'southwest', 'north', 'south', 'east', 'west'};
            app.width_box = 1.8 * app.width_box;
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, list_legendlocation)
            app.width_box = app.width_box / 1.8;
        end

        function fill_plotConfig_axes_node(app)
            % Fill plotConfigs_axes_node

            % Initialization
            obj_num = 0;
            obj_num_children = 0;

            % Title
            create_title(app, 'Combustion Toolbox plot configuration: Axes')

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
            tag = {'plotConfig', 'box'};
            title = 'Display axes outline';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'})

            % ** Miscellaneous: Axes - grid
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'grid'};
            title = 'Display or hide axes grid lines';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'})

            % ** Miscellaneous: Axes - hold
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'hold'};
            title = 'Retain current plot when adding new plots';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'on', 'off'})

            % ** Miscellaneous: Axes - axis_x
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'axis_x'};
            title = 'Set x-axis limits';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'auto', 'tight', 'padded', 'fill', 'equal', 'image', 'square', 'vis3d'})

            % ** Miscellaneous: Axes - axis_y
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'axis_y'};
            title = 'Set y-axis limits';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'auto', 'tight', 'padded', 'fill', 'equal', 'image', 'square', 'vis3d'})

            % ** Miscellaneous: Axes - xscale
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'xscale'};
            title = 'Set x-axis scale';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'linear', 'log'})

            % ** Miscellaneous: Axes - yscale
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'yscale'};
            title = 'Set y-axis scale';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'linear', 'log'})

            % ** Miscellaneous: Axes - xdir
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'xdir'};
            title = 'Set x-axis direction';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'normal', 'reverse'})

            % ** Miscellaneous: Axes - ydir
            obj_num_children = obj_num_children + 1;
            tag = {'plotConfig', 'ydir'};
            title = 'Set y-axis direction';
            create_panel_dropdown(app, panel_name, obj_num_children, tag, title, {'normal', 'reverse'})
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
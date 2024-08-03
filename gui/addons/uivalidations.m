classdef uivalidations < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIValidations    matlab.ui.Figure
        text_phi_4       matlab.ui.control.Label
        matlab           matlab.ui.control.EditField
        Lamp             matlab.ui.control.Lamp
        problems_solved  matlab.ui.control.NumericEditField
        CPUtimesLabel_4  matlab.ui.control.Label
        text_phi_3       matlab.ui.control.Label
        OS               matlab.ui.control.EditField
        cores            matlab.ui.control.NumericEditField
        CPUtimesLabel_3  matlab.ui.control.Label
        memory           matlab.ui.control.NumericEditField
        CPUtimesLabel_2  matlab.ui.control.Label
        text_phi         matlab.ui.control.Label
        cpu              matlab.ui.control.EditField
        cputime          matlab.ui.control.NumericEditField
        CPUtimesLabel    matlab.ui.control.Label
        CloseButton      matlab.ui.control.Button
        RunButton        matlab.ui.control.Button
        Tree             matlab.ui.container.Tree
        CEA              matlab.ui.container.TreeNode
        SDToolbox        matlab.ui.container.TreeNode
        CANTERA          matlab.ui.container.TreeNode
        TEA              matlab.ui.container.TreeNode
    end

    
    properties (Access = private)
        filename = 'LOG.txt'; % Filename of the log file
        system                % Info of the system
        color_lamp_working = [0.9961, 0.9804, 0.8314]; % Lamp color (rgb): working
        color_lamp_done    = [0.5608, 0.7255, 0.6588]; % Lamp color (rgb): done
        color_lamp_error   = [0.9451, 0.5059, 0.5529]; % Lamp color (rgb): error
    end
    
    methods (Access = private)
        
        function start_log(app)
            if exist(app.filename, 'file') 
                delete(app.filename); 
            end
            diary(app.filename);
            diary on
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Load info system
            app.system = combustiontoolbox.utils.extensions.cpuinfo();
            % Assign info of the system
            app.OS.Value = [app.system.OSType, ' ', app.system.OSVersion];
            app.cpu.Value = app.system.CPUName;
            app.matlab.Value = ['Release R' version('-release')];
            app.memory.Value = app.system.TotalMemory * 9.313225746154785e-10; % [GB]
            app.cores.Value = app.system.TotalCores;
            % Include validations in the UITree object
            gui_add_nodes_validations(app, 'CEA');
            gui_add_nodes_validations(app, 'SDToolbox');
            gui_add_nodes_validations(app, 'CANTERA');
            gui_add_nodes_validations(app, 'TEA');
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            try
                % Initialization
                close all;
                app.cputime.Value = 0;
                app.problems_solved.Value = 0;
                app.Lamp.Color = app.color_lamp_working;
                % CPU time
                cputime = tic;
                % Run validation
                % app.problems_solved.Value = eval(app.Tree.SelectedNodes.Text(1:end-2));
                eval(app.Tree.SelectedNodes.Text(1:end-2));
                % CPU time
                app.cputime.Value = toc(cputime);
                % Done!
                app.Lamp.Color = app.color_lamp_done;
            catch ME
                app.Lamp.Color = app.color_lamp_error;
                % Print error
                message = {sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                ME.stack(1).name, ME.stack(1).line, ME.message)};
                uialert(app.UIValidations, message, 'Warning', 'Icon', 'warning');
            end
        end

        % Close request function: UIValidations
        function UIValidationsCloseRequest(app, event)
            delete(app)
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            UIValidationsCloseRequest(app, event)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIValidations and hide until all components are created
            app.UIValidations = uifigure('Visible', 'off');
            app.UIValidations.AutoResizeChildren = 'off';
            app.UIValidations.Color = [0.9098 0.9098 0.8902];
            app.UIValidations.Position = [650 400 670 281];
            app.UIValidations.Name = 'UIValidations';
            app.UIValidations.Icon = 'C:\Users\Alber\Mi unidad\Phd\Combustion_Toolbox\Combustion_Toolbox_repo\combustion_toolbox\gui\assets\uivalidations\uivalidation.png';
            app.UIValidations.Resize = 'off';
            app.UIValidations.CloseRequestFcn = createCallbackFcn(app, @UIValidationsCloseRequest, true);
            app.UIValidations.Scrollable = 'on';

            % Create Tree
            app.Tree = uitree(app.UIValidations);
            app.Tree.Position = [14 15 308 253];

            % Create CEA
            app.CEA = uitreenode(app.Tree);
            app.CEA.Text = 'CEA';

            % Create SDToolbox
            app.SDToolbox = uitreenode(app.Tree);
            app.SDToolbox.Text = 'SDToolbox';

            % Create CANTERA
            app.CANTERA = uitreenode(app.Tree);
            app.CANTERA.Text = 'CANTERA';

            % Create TEA
            app.TEA = uitreenode(app.Tree);
            app.TEA.Text = 'TEA';

            % Create RunButton
            app.RunButton = uibutton(app.UIValidations, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Position = [417 15 90 22];
            app.RunButton.Text = 'Run';

            % Create CloseButton
            app.CloseButton = uibutton(app.UIValidations, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Position = [525 15 90 22];
            app.CloseButton.Text = 'Close';

            % Create CPUtimesLabel
            app.CPUtimesLabel = uilabel(app.UIValidations);
            app.CPUtimesLabel.HorizontalAlignment = 'right';
            app.CPUtimesLabel.Position = [429 88 69 22];
            app.CPUtimesLabel.Text = 'CPU time';

            % Create cputime
            app.cputime = uieditfield(app.UIValidations, 'numeric');
            app.cputime.ValueDisplayFormat = '%11.4g s';
            app.cputime.Editable = 'off';
            app.cputime.Position = [508 88 147 22];

            % Create cpu
            app.cpu = uieditfield(app.UIValidations, 'text');
            app.cpu.Editable = 'off';
            app.cpu.HorizontalAlignment = 'center';
            app.cpu.Position = [370 212 285 22];

            % Create text_phi
            app.text_phi = uilabel(app.UIValidations);
            app.text_phi.HorizontalAlignment = 'right';
            app.text_phi.Position = [321 217 39 25];
            app.text_phi.Text = 'CPU';

            % Create CPUtimesLabel_2
            app.CPUtimesLabel_2 = uilabel(app.UIValidations);
            app.CPUtimesLabel_2.HorizontalAlignment = 'right';
            app.CPUtimesLabel_2.Position = [420 148 78 22];
            app.CPUtimesLabel_2.Text = 'Total Memory';

            % Create memory
            app.memory = uieditfield(app.UIValidations, 'numeric');
            app.memory.ValueDisplayFormat = '%11.4g GB';
            app.memory.Editable = 'off';
            app.memory.Position = [508 148 147 22];

            % Create CPUtimesLabel_3
            app.CPUtimesLabel_3 = uilabel(app.UIValidations);
            app.CPUtimesLabel_3.HorizontalAlignment = 'right';
            app.CPUtimesLabel_3.Position = [460 118 38 22];
            app.CPUtimesLabel_3.Text = 'Cores';

            % Create cores
            app.cores = uieditfield(app.UIValidations, 'numeric');
            app.cores.Editable = 'off';
            app.cores.Position = [508 118 147 22];

            % Create OS
            app.OS = uieditfield(app.UIValidations, 'text');
            app.OS.Editable = 'off';
            app.OS.HorizontalAlignment = 'center';
            app.OS.Position = [370 245 285 22];

            % Create text_phi_3
            app.text_phi_3 = uilabel(app.UIValidations);
            app.text_phi_3.HorizontalAlignment = 'right';
            app.text_phi_3.Position = [321 243 39 25];
            app.text_phi_3.Text = 'OS';

            % Create CPUtimesLabel_4
            app.CPUtimesLabel_4 = uilabel(app.UIValidations);
            app.CPUtimesLabel_4.HorizontalAlignment = 'right';
            app.CPUtimesLabel_4.Enable = 'off';
            app.CPUtimesLabel_4.Visible = 'off';
            app.CPUtimesLabel_4.Position = [404 58 94 22];
            app.CPUtimesLabel_4.Text = 'Problems solved';

            % Create problems_solved
            app.problems_solved = uieditfield(app.UIValidations, 'numeric');
            app.problems_solved.Editable = 'off';
            app.problems_solved.Enable = 'off';
            app.problems_solved.Visible = 'off';
            app.problems_solved.Position = [508 58 147 22];

            % Create Lamp
            app.Lamp = uilamp(app.UIValidations);
            app.Lamp.Position = [632 15 23 23];
            app.Lamp.Color = [0.8 0.8 0.8];

            % Create matlab
            app.matlab = uieditfield(app.UIValidations, 'text');
            app.matlab.Editable = 'off';
            app.matlab.HorizontalAlignment = 'center';
            app.matlab.Position = [508 180 147 22];

            % Create text_phi_4
            app.text_phi_4 = uilabel(app.UIValidations);
            app.text_phi_4.HorizontalAlignment = 'right';
            app.text_phi_4.Position = [417 178 81 25];
            app.text_phi_4.Text = 'MATLAB';

            % Show the figure after all components are created
            app.UIValidations.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = uivalidations

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIValidations)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UIValidations)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIValidations)
        end
    end
end
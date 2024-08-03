classdef uilicense < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UILicense    matlab.ui.Figure
        Hyperlink    matlab.ui.control.Hyperlink
        Label_3      matlab.ui.control.Label
        Label_2      matlab.ui.control.Label
        Label_title  matlab.ui.control.Label
        OkButton     matlab.ui.control.Button
        TextArea     matlab.ui.control.TextArea
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Set version Combustion Toolbox
            release = combustiontoolbox.common.Constants.release;
            app.Label_title.Text = ['Combustion Toolbox ', release];
            % Set license
            app.TextArea.Value = GPL();
        end

        % Button pushed function: OkButton
        function OkButtonPushed(app, event)
            UILicenseCloseRequest(app, event);
        end

        % Close request function: UILicense
        function UILicenseCloseRequest(app, event)
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UILicense and hide until all components are created
            app.UILicense = uifigure('Visible', 'off');
            app.UILicense.Color = [0.9098 0.9098 0.8902];
            app.UILicense.Position = [650 250 618 439];
            app.UILicense.Name = 'License - Combustion Toolbox';
            app.UILicense.Icon = 'logo_uilicense.jpg';
            app.UILicense.CloseRequestFcn = createCallbackFcn(app, @UILicenseCloseRequest, true);

            % Create TextArea
            app.TextArea = uitextarea(app.UILicense);
            app.TextArea.Position = [16 42 588 217];

            % Create OkButton
            app.OkButton = uibutton(app.UILicense, 'push');
            app.OkButton.ButtonPushedFcn = createCallbackFcn(app, @OkButtonPushed, true);
            app.OkButton.Position = [504 9 100 23];
            app.OkButton.Text = 'Ok';

            % Create Label_title
            app.Label_title = uilabel(app.UILicense);
            app.Label_title.HorizontalAlignment = 'center';
            app.Label_title.FontSize = 14;
            app.Label_title.FontWeight = 'bold';
            app.Label_title.FontColor = [0.6353 0.0784 0.1843];
            app.Label_title.Position = [111 408 394 22];
            app.Label_title.Text = 'Combustion Toolbox vXX.XX.XX';

            % Create Label_2
            app.Label_2 = uilabel(app.UILicense);
            app.Label_2.HorizontalAlignment = 'center';
            app.Label_2.Position = [111 376 399 30];
            app.Label_2.Text = {'A MATLAB-GUI based open-source tool for solving combustion problems'; sprintf('Copyright (C) 2022-%d Alberto Cuadra Lara and contributors', year(datetime('now')))};

            % Create Label_3
            app.Label_3 = uilabel(app.UILicense);
            app.Label_3.HorizontalAlignment = 'center';
            app.Label_3.Position = [101 273 419 59];
            app.Label_3.Text = {'This program is distributed in the hope that it will be useful,'; 'but WITHOUT ANY WARRANTY; without even the implied warranty of'; 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the'; 'GNU General Public License for more details.'};

            % Create Hyperlink
            app.Hyperlink = uihyperlink(app.UILicense);
            app.Hyperlink.HorizontalAlignment = 'center';
            app.Hyperlink.URL = combustiontoolbox.utils.SystemUtils.url.website;
            app.Hyperlink.Position = [158 355 304 22];
            app.Hyperlink.Text = combustiontoolbox.utils.SystemUtils.url.website;

            % Show the figure after all components are created
            app.UILicense.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = uilicense

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UILicense)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UILicense)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UILicense)
        end
    end
end
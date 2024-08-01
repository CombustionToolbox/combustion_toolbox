classdef uiabout < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIAbout            matlab.ui.Figure
        Panel_bar          matlab.ui.container.Panel
        GridLayout         matlab.ui.container.GridLayout
        Image_logo_1       matlab.ui.control.Image
        Image_logo_2       matlab.ui.control.Image
        Image_bar_1        matlab.ui.control.Image
        Image_bar_2        matlab.ui.control.Image
        Image_bar_3        matlab.ui.control.Image
        Panel_authors      matlab.ui.container.Panel
        GridLayout2        matlab.ui.container.GridLayout
        Panel_acuadra      matlab.ui.container.Panel
        GridLayoutacuadra  matlab.ui.container.GridLayout
        acuadraImage       matlab.ui.control.Image
        acuadraMail        matlab.ui.control.Image
        acuadraLinkedin    matlab.ui.control.Image
        acuadraResearch    matlab.ui.control.Image
        acuadraGithub      matlab.ui.control.Image
        categLabel1        matlab.ui.control.Label
        nameLabel1         matlab.ui.control.Label
        Panel_mvera        matlab.ui.container.Panel
        GridLayout_mvera   matlab.ui.container.GridLayout
        mveraImage         matlab.ui.control.Image
        mveraLinkedin      matlab.ui.control.Image
        mveraResearch      matlab.ui.control.Image
        nameLabel2         matlab.ui.control.Label
        categLabel2        matlab.ui.control.Label
        Panel_chuete       matlab.ui.container.Panel
        GridLayout_chuete  matlab.ui.container.GridLayout
        chueteImage        matlab.ui.control.Image
        chueteLinkedin     matlab.ui.control.Image
        chueteResearch     matlab.ui.control.Image
        categLabel3        matlab.ui.control.Label
        nameLabel3         matlab.ui.control.Label
    end

    
    methods (Access = private)
        
        function select_action(app, event)
            % Import packages
            import combustiontoolbox.common.Constants
            import combustiontoolbox.utils.openWebsite

            obj_tag_clicked = event.Source.Tag;
            switch lower(obj_tag_clicked)
                case {'repo', 'repository'}
                    property = 'repository';
                case 'website_ct'
                    property = 'website';
                case {'docs', 'documentation'}
                    property = 'documentation';
                case 'publications'
                    property = 'publications';
                case 'website_uc3m'
                    property = 'websiteUC3M';
                case 'website_acuadra'
                    property = 'websiteACuadra';
                case 'mail_acuadra'
                    property = 'mailACuadra';
                case 'github_acuadra'
                    property = 'githubACuadra';
                case 'researchgate_acuadra'
                    property = 'researchgateACuadra';
                case 'linkedin_acuadra'
                    property = 'linkedinACuadra';
                case 'website_mvera'
                    property = 'websiteMVera';
                case 'researchgate_mvera'
                    property = 'researchgateMVera';
                case 'linkedin_mvera'
                    property = 'linkedinMVera';
                case 'website_chuete'
                    property = 'websiteCHuete';
                case 'researchgate_chuete'
                    property = 'researchgateCHuete';
                case 'linkedin_chuete'
                    property = 'linkedinCHuete';
                case 'contributors'
                    property = 'contributors';
            end
            
            url = Constants.url.(property);
            openWebsite(url);
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
%             gui_SnapshotMenuSelected(app.UIAbout);
        end

        % Image clicked function: Image_bar_1, Image_bar_2, Image_bar_3, 
        % ...and 13 other components
        function Image_2Clicked(app, event)
            select_action(app, event);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIAbout and hide until all components are created
            app.UIAbout = uifigure('Visible', 'off');
            app.UIAbout.Color = [0.9098 0.9098 0.8902];
            app.UIAbout.Position = [650 300 699 350];
            app.UIAbout.Name = 'About';
            app.UIAbout.Icon = 'logo_uiabout.png';
            app.UIAbout.Scrollable = 'on';

            % Create Panel_authors
            app.Panel_authors = uipanel(app.UIAbout);
            app.Panel_authors.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_authors.Position = [0 0 700 350];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Panel_authors);
            app.GridLayout2.ColumnWidth = {'1x', '1x', '1x'};
            app.GridLayout2.RowHeight = {36, '1x'};
            app.GridLayout2.ColumnSpacing = 0;
            app.GridLayout2.Padding = [0 0 0 0];
            app.GridLayout2.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create Panel_chuete
            app.Panel_chuete = uipanel(app.GridLayout2);
            app.Panel_chuete.BorderType = 'none';
            app.Panel_chuete.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_chuete.Layout.Row = 2;
            app.Panel_chuete.Layout.Column = 3;

            % Create GridLayout_chuete
            app.GridLayout_chuete = uigridlayout(app.Panel_chuete);
            app.GridLayout_chuete.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout_chuete.RowHeight = {'5x', '1x', '1x', '1.1x'};
            app.GridLayout_chuete.ColumnSpacing = 0;
            app.GridLayout_chuete.RowSpacing = 0;
            app.GridLayout_chuete.Padding = [5 5 5 5];
            app.GridLayout_chuete.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create nameLabel3
            app.nameLabel3 = uilabel(app.GridLayout_chuete);
            app.nameLabel3.HorizontalAlignment = 'center';
            app.nameLabel3.VerticalAlignment = 'bottom';
            app.nameLabel3.FontSize = 18;
            app.nameLabel3.FontWeight = 'bold';
            app.nameLabel3.FontColor = [0.3412 0.4784 0.5843];
            app.nameLabel3.Layout.Row = 2;
            app.nameLabel3.Layout.Column = [1 4];
            app.nameLabel3.Text = 'César Huete';

            % Create categLabel3
            app.categLabel3 = uilabel(app.GridLayout_chuete);
            app.categLabel3.HorizontalAlignment = 'center';
            app.categLabel3.VerticalAlignment = 'top';
            app.categLabel3.FontSize = 16;
            app.categLabel3.FontColor = [0.6 0.6 0.6];
            app.categLabel3.Layout.Row = 3;
            app.categLabel3.Layout.Column = [1 4];
            app.categLabel3.Text = 'Advisor';

            % Create chueteResearch
            app.chueteResearch = uiimage(app.GridLayout_chuete);
            app.chueteResearch.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.chueteResearch.Tag = 'researchgate_chuete';
            app.chueteResearch.Layout.Row = 4;
            app.chueteResearch.Layout.Column = 2;
            app.chueteResearch.ImageSource = 'icon_researchgate.svg';

            % Create chueteLinkedin
            app.chueteLinkedin = uiimage(app.GridLayout_chuete);
            app.chueteLinkedin.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.chueteLinkedin.Tag = 'linkedin_chuete';
            app.chueteLinkedin.Layout.Row = 4;
            app.chueteLinkedin.Layout.Column = 3;
            app.chueteLinkedin.ImageSource = 'icon_linkedin.svg';

            % Create chueteImage
            app.chueteImage = uiimage(app.GridLayout_chuete);
            app.chueteImage.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.chueteImage.Tag = 'website_chuete';
            app.chueteImage.Layout.Row = 1;
            app.chueteImage.Layout.Column = [1 4];
            app.chueteImage.ImageSource = 'profile_chuete.svg';

            % Create Panel_mvera
            app.Panel_mvera = uipanel(app.GridLayout2);
            app.Panel_mvera.BorderType = 'none';
            app.Panel_mvera.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_mvera.Layout.Row = 2;
            app.Panel_mvera.Layout.Column = 2;

            % Create GridLayout_mvera
            app.GridLayout_mvera = uigridlayout(app.Panel_mvera);
            app.GridLayout_mvera.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayout_mvera.RowHeight = {'5x', '1x', '1x', '1.1x'};
            app.GridLayout_mvera.ColumnSpacing = 0;
            app.GridLayout_mvera.RowSpacing = 0;
            app.GridLayout_mvera.Padding = [5 5 5 5];
            app.GridLayout_mvera.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create categLabel2
            app.categLabel2 = uilabel(app.GridLayout_mvera);
            app.categLabel2.HorizontalAlignment = 'center';
            app.categLabel2.VerticalAlignment = 'top';
            app.categLabel2.FontSize = 16;
            app.categLabel2.FontColor = [0.6 0.6 0.6];
            app.categLabel2.Layout.Row = 3;
            app.categLabel2.Layout.Column = [1 4];
            app.categLabel2.Text = 'Advisor';

            % Create nameLabel2
            app.nameLabel2 = uilabel(app.GridLayout_mvera);
            app.nameLabel2.HorizontalAlignment = 'center';
            app.nameLabel2.VerticalAlignment = 'bottom';
            app.nameLabel2.FontSize = 18;
            app.nameLabel2.FontWeight = 'bold';
            app.nameLabel2.FontColor = [0.3412 0.4784 0.5843];
            app.nameLabel2.Layout.Row = 2;
            app.nameLabel2.Layout.Column = [1 4];
            app.nameLabel2.Text = 'Marcos Vera';

            % Create mveraResearch
            app.mveraResearch = uiimage(app.GridLayout_mvera);
            app.mveraResearch.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.mveraResearch.Tag = 'researchgate_mvera';
            app.mveraResearch.Layout.Row = 4;
            app.mveraResearch.Layout.Column = 2;
            app.mveraResearch.ImageSource = 'icon_researchgate.svg';

            % Create mveraLinkedin
            app.mveraLinkedin = uiimage(app.GridLayout_mvera);
            app.mveraLinkedin.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.mveraLinkedin.Tag = 'linkedin_mvera';
            app.mveraLinkedin.Layout.Row = 4;
            app.mveraLinkedin.Layout.Column = 3;
            app.mveraLinkedin.ImageSource = 'icon_linkedin.svg';

            % Create mveraImage
            app.mveraImage = uiimage(app.GridLayout_mvera);
            app.mveraImage.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.mveraImage.Tag = 'website_mvera';
            app.mveraImage.Layout.Row = 1;
            app.mveraImage.Layout.Column = [1 4];
            app.mveraImage.ImageSource = 'profile_mvera.svg';

            % Create Panel_acuadra
            app.Panel_acuadra = uipanel(app.GridLayout2);
            app.Panel_acuadra.BorderType = 'none';
            app.Panel_acuadra.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_acuadra.Layout.Row = 2;
            app.Panel_acuadra.Layout.Column = 1;

            % Create GridLayoutacuadra
            app.GridLayoutacuadra = uigridlayout(app.Panel_acuadra);
            app.GridLayoutacuadra.ColumnWidth = {'1x', '1x', '1x', '1x'};
            app.GridLayoutacuadra.RowHeight = {'5x', '1x', '1x', '1.1x'};
            app.GridLayoutacuadra.ColumnSpacing = 0;
            app.GridLayoutacuadra.RowSpacing = 0;
            app.GridLayoutacuadra.Padding = [5 5 5 5];
            app.GridLayoutacuadra.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create nameLabel1
            app.nameLabel1 = uilabel(app.GridLayoutacuadra);
            app.nameLabel1.HorizontalAlignment = 'center';
            app.nameLabel1.VerticalAlignment = 'bottom';
            app.nameLabel1.FontSize = 18;
            app.nameLabel1.FontWeight = 'bold';
            app.nameLabel1.FontColor = [0.3412 0.4784 0.5804];
            app.nameLabel1.Layout.Row = 2;
            app.nameLabel1.Layout.Column = [1 4];
            app.nameLabel1.Text = 'Alberto Cuadra-Lara';

            % Create categLabel1
            app.categLabel1 = uilabel(app.GridLayoutacuadra);
            app.categLabel1.HorizontalAlignment = 'center';
            app.categLabel1.VerticalAlignment = 'top';
            app.categLabel1.FontSize = 16;
            app.categLabel1.FontColor = [0.6 0.6 0.6];
            app.categLabel1.Layout.Row = 3;
            app.categLabel1.Layout.Column = [1 4];
            app.categLabel1.Text = 'Lead developer';

            % Create acuadraGithub
            app.acuadraGithub = uiimage(app.GridLayoutacuadra);
            app.acuadraGithub.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.acuadraGithub.Tag = 'github_acuadra';
            app.acuadraGithub.Layout.Row = 4;
            app.acuadraGithub.Layout.Column = 2;
            app.acuadraGithub.ImageSource = 'icon_github.svg';

            % Create acuadraResearch
            app.acuadraResearch = uiimage(app.GridLayoutacuadra);
            app.acuadraResearch.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.acuadraResearch.Tag = 'researchgate_acuadra';
            app.acuadraResearch.Layout.Row = 4;
            app.acuadraResearch.Layout.Column = 3;
            app.acuadraResearch.ImageSource = 'icon_researchgate.svg';

            % Create acuadraLinkedin
            app.acuadraLinkedin = uiimage(app.GridLayoutacuadra);
            app.acuadraLinkedin.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.acuadraLinkedin.Tag = 'linkedin_acuadra';
            app.acuadraLinkedin.Layout.Row = 4;
            app.acuadraLinkedin.Layout.Column = 4;
            app.acuadraLinkedin.ImageSource = 'icon_linkedin.svg';

            % Create acuadraMail
            app.acuadraMail = uiimage(app.GridLayoutacuadra);
            app.acuadraMail.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.acuadraMail.Tag = 'mail_acuadra';
            app.acuadraMail.Layout.Row = 4;
            app.acuadraMail.Layout.Column = 1;
            app.acuadraMail.ImageSource = 'icon_mail.svg';

            % Create acuadraImage
            app.acuadraImage = uiimage(app.GridLayoutacuadra);
            app.acuadraImage.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.acuadraImage.Tag = 'website_acuadra';
            app.acuadraImage.Layout.Row = 1;
            app.acuadraImage.Layout.Column = [1 4];
            app.acuadraImage.ImageSource = 'profile_acuadra.svg';

            % Create Panel_bar
            app.Panel_bar = uipanel(app.UIAbout);
            app.Panel_bar.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_bar.Position = [0 307 700 44];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Panel_bar);
            app.GridLayout.ColumnWidth = {36, 142, 99, 76, '1x', 190};
            app.GridLayout.RowHeight = {36};
            app.GridLayout.RowSpacing = 3.2400016784668;
            app.GridLayout.Padding = [5 3.2400016784668 5 3.2400016784668];
            app.GridLayout.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create Image_bar_3
            app.Image_bar_3 = uiimage(app.GridLayout);
            app.Image_bar_3.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.Image_bar_3.Tag = 'publications';
            app.Image_bar_3.Layout.Row = 1;
            app.Image_bar_3.Layout.Column = 4;
            app.Image_bar_3.ImageSource = 'bar_publications.svg';

            % Create Image_bar_2
            app.Image_bar_2 = uiimage(app.GridLayout);
            app.Image_bar_2.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.Image_bar_2.Tag = 'documentation';
            app.Image_bar_2.Layout.Row = 1;
            app.Image_bar_2.Layout.Column = 3;
            app.Image_bar_2.ImageSource = 'bar_documentation.svg';

            % Create Image_bar_1
            app.Image_bar_1 = uiimage(app.GridLayout);
            app.Image_bar_1.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.Image_bar_1.Tag = 'website_CT';
            app.Image_bar_1.Layout.Row = 1;
            app.Image_bar_1.Layout.Column = 2;
            app.Image_bar_1.ImageSource = 'bar_home.svg';

            % Create Image_logo_2
            app.Image_logo_2 = uiimage(app.GridLayout);
            app.Image_logo_2.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.Image_logo_2.Tag = 'website_uc3m';
            app.Image_logo_2.Layout.Row = 1;
            app.Image_logo_2.Layout.Column = 6;
            app.Image_logo_2.HorizontalAlignment = 'right';
            app.Image_logo_2.ImageSource = 'logo_fluids_uc3m.svg';

            % Create Image_logo_1
            app.Image_logo_1 = uiimage(app.GridLayout);
            app.Image_logo_1.ImageClickedFcn = createCallbackFcn(app, @Image_2Clicked, true);
            app.Image_logo_1.Tag = 'website_CT';
            app.Image_logo_1.Layout.Row = 1;
            app.Image_logo_1.Layout.Column = 1;
            app.Image_logo_1.ImageSource = 'logo_CT_noversion.svg';

            % Show the figure after all components are created
            app.UIAbout.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = uiabout

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIAbout)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UIAbout)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIAbout)
        end
    end
end
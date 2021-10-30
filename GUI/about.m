classdef about < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AboutUIFigure       matlab.ui.Figure
        fluidosImage        matlab.ui.control.Image
        chueteImage         matlab.ui.control.Image
        mveraImage          matlab.ui.control.Image
        acuadraImage        matlab.ui.control.Image
        chueteLinkedin      matlab.ui.control.Image
        chueteResearch      matlab.ui.control.Image
        mveraLinkedin       matlab.ui.control.Image
        mveraResearch       matlab.ui.control.Image
        acuadraMail         matlab.ui.control.Image
        acuadraLinkedin     matlab.ui.control.Image
        acuadraResearch     matlab.ui.control.Image
        acuadraGithub       matlab.ui.control.Image
        categLabel3         matlab.ui.control.Label
        nameLabel3          matlab.ui.control.Label
        nameLabel2          matlab.ui.control.Label
        categLabel2         matlab.ui.control.Label
        RepositoryButton    matlab.ui.control.Button
        ContributorsButton  matlab.ui.control.Button
        WikiButton          matlab.ui.control.Button
        categLabel1         matlab.ui.control.Label
        nameLabel1          matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
%             gui_SnapshotMenuSelected(app.AboutUIFigure);
        end

        % Image clicked function: fluidosImage
        function fluidosWebClicked(app, event)
            web('http://fluidosuc3m.es/','-browser');
        end

        % Image clicked function: acuadraLinkedin
        function acuadraLinkedinClicked(app, event)
            web('https://www.linkedin.com/in/albertocuadralara/','-browser');
        end

        % Image clicked function: mveraLinkedin
        function mveraLinkedinClicked(app, event)
            web('https://www.linkedin.com/in/marcos-vera-coello-05b3a643/?originalSubdomain=es','-browser');
        end

        % Image clicked function: acuadraMail
        function acuadraMailClicked(app, event)
            web('mailto:acuadra@ing.uc3m.es');
        end

        % Callback function
        function Button_9Pushed(app, event)
            
        end

        % Image clicked function: mveraResearch
        function mveraResearchClicked(app, event)
             web('https://www.researchgate.net/profile/Marcos_Vera','-browser');
        end

        % Image clicked function: acuadraResearch
        function acuadraResearchClicked(app, event)
            web('https://www.researchgate.net/profile/Alberto_Cuadra_Lara','-browser');
        end

        % Image clicked function: acuadraImage
        function acuadraWebClicked(app, event)
            web('https://acuadralara.com/','-browser');
        end

        % Image clicked function: mveraImage
        function mveraWebClicked(app, event)
            web('http://fluidosuc3m.es/people/mvcoello/','-browser');
        end

        % Image clicked function: acuadraGithub
        function acuadraGithubClicked(app, event)
            web('https://github.com/AlbertoCuadra','-browser');
        end

        % Button pushed function: RepositoryButton
        function RepositoryButtonPushed(app, event)
            web('https://github.com/AlbertoCuadra/combustion_toolbox','-browser');
        end

        % Button pushed function: ContributorsButton
        function ContributorsButtonPushed(app, event)
            web('https://github.com/AlbertoCuadra/combustion_toolbox/blob/master/CONTRIBUTORS.md','-browser');
        end

        % Image clicked function: chueteResearch
        function chueteResearchClicked(app, event)
            web('https://www.researchgate.net/profile/Cesar-Huete','-browser');
        end

        % Image clicked function: chueteLinkedin
        function chueteLinkedinClicked(app, event)
            web('https://www.linkedin.com/in/cesarhuete/','-browser');
        end

        % Image clicked function: chueteImage
        function chueteWebClicked(app, event)
            web('http://fluidosuc3m.es/people/chuete/','-browser');
        end

        % Button pushed function: WikiButton
        function WikiButtonPushed(app, event)
            web('https://github.com/AlbertoCuadra/combustion_toolbox/wiki','-browser');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AboutUIFigure and hide until all components are created
            app.AboutUIFigure = uifigure('Visible', 'off');
            app.AboutUIFigure.AutoResizeChildren = 'off';
            app.AboutUIFigure.Color = [0.9098 0.9098 0.8941];
            app.AboutUIFigure.Position = [650 300 634 389];
            app.AboutUIFigure.Name = 'About';
            app.AboutUIFigure.Resize = 'off';

            % Create nameLabel1
            app.nameLabel1 = uilabel(app.AboutUIFigure);
            app.nameLabel1.HorizontalAlignment = 'center';
            app.nameLabel1.FontSize = 18;
            app.nameLabel1.FontWeight = 'bold';
            app.nameLabel1.FontColor = [0.3412 0.4784 0.5804];
            app.nameLabel1.Position = [23 188 207 22];
            app.nameLabel1.Text = 'Alberto Cuadra-Lara';

            % Create categLabel1
            app.categLabel1 = uilabel(app.AboutUIFigure);
            app.categLabel1.HorizontalAlignment = 'center';
            app.categLabel1.FontSize = 16;
            app.categLabel1.FontColor = [0.6 0.6 0.6];
            app.categLabel1.Position = [8 161 236 22];
            app.categLabel1.Text = 'Core developer & App Designer';

            % Create WikiButton
            app.WikiButton = uibutton(app.AboutUIFigure, 'push');
            app.WikiButton.ButtonPushedFcn = createCallbackFcn(app, @WikiButtonPushed, true);
            app.WikiButton.Position = [495 14 127 22];
            app.WikiButton.Text = 'Wiki';

            % Create ContributorsButton
            app.ContributorsButton = uibutton(app.AboutUIFigure, 'push');
            app.ContributorsButton.ButtonPushedFcn = createCallbackFcn(app, @ContributorsButtonPushed, true);
            app.ContributorsButton.Position = [495 41 127 22];
            app.ContributorsButton.Text = 'Contributors';

            % Create RepositoryButton
            app.RepositoryButton = uibutton(app.AboutUIFigure, 'push');
            app.RepositoryButton.ButtonPushedFcn = createCallbackFcn(app, @RepositoryButtonPushed, true);
            app.RepositoryButton.Position = [495 68 127 22];
            app.RepositoryButton.Text = 'Repository';

            % Create categLabel2
            app.categLabel2 = uilabel(app.AboutUIFigure);
            app.categLabel2.HorizontalAlignment = 'center';
            app.categLabel2.FontSize = 16;
            app.categLabel2.FontColor = [0.6 0.6 0.6];
            app.categLabel2.Position = [283 161 109 22];
            app.categLabel2.Text = 'Developer';

            % Create nameLabel2
            app.nameLabel2 = uilabel(app.AboutUIFigure);
            app.nameLabel2.HorizontalAlignment = 'center';
            app.nameLabel2.FontSize = 18;
            app.nameLabel2.FontWeight = 'bold';
            app.nameLabel2.FontColor = [0.3412 0.4784 0.5843];
            app.nameLabel2.Position = [283 188 109 22];
            app.nameLabel2.Text = 'Marcos Vera';

            % Create nameLabel3
            app.nameLabel3 = uilabel(app.AboutUIFigure);
            app.nameLabel3.HorizontalAlignment = 'center';
            app.nameLabel3.FontSize = 18;
            app.nameLabel3.FontWeight = 'bold';
            app.nameLabel3.FontColor = [0.3412 0.4784 0.5843];
            app.nameLabel3.Position = [488 188 127 22];
            app.nameLabel3.Text = 'César Huete';

            % Create categLabel3
            app.categLabel3 = uilabel(app.AboutUIFigure);
            app.categLabel3.HorizontalAlignment = 'center';
            app.categLabel3.FontSize = 16;
            app.categLabel3.FontColor = [0.6 0.6 0.6];
            app.categLabel3.Position = [488 161 127 22];
            app.categLabel3.Text = 'Developer';

            % Create acuadraGithub
            app.acuadraGithub = uiimage(app.AboutUIFigure);
            app.acuadraGithub.ImageClickedFcn = createCallbackFcn(app, @acuadraGithubClicked, true);
            app.acuadraGithub.Position = [84 120 43 37];
            app.acuadraGithub.ImageSource = 'github_icon.svg';

            % Create acuadraResearch
            app.acuadraResearch = uiimage(app.AboutUIFigure);
            app.acuadraResearch.ImageClickedFcn = createCallbackFcn(app, @acuadraResearchClicked, true);
            app.acuadraResearch.Position = [133 120 43 37];
            app.acuadraResearch.ImageSource = 'researchgate_icon.svg';

            % Create acuadraLinkedin
            app.acuadraLinkedin = uiimage(app.AboutUIFigure);
            app.acuadraLinkedin.ImageClickedFcn = createCallbackFcn(app, @acuadraLinkedinClicked, true);
            app.acuadraLinkedin.Position = [182 120 43 37];
            app.acuadraLinkedin.ImageSource = 'linkedin_icon.svg';

            % Create acuadraMail
            app.acuadraMail = uiimage(app.AboutUIFigure);
            app.acuadraMail.ImageClickedFcn = createCallbackFcn(app, @acuadraMailClicked, true);
            app.acuadraMail.Position = [35 120 43 37];
            app.acuadraMail.ImageSource = 'mail_icon.svg';

            % Create mveraResearch
            app.mveraResearch = uiimage(app.AboutUIFigure);
            app.mveraResearch.ImageClickedFcn = createCallbackFcn(app, @mveraResearchClicked, true);
            app.mveraResearch.Position = [293 120 43 37];
            app.mveraResearch.ImageSource = 'researchgate_icon.svg';

            % Create mveraLinkedin
            app.mveraLinkedin = uiimage(app.AboutUIFigure);
            app.mveraLinkedin.ImageClickedFcn = createCallbackFcn(app, @mveraLinkedinClicked, true);
            app.mveraLinkedin.Position = [342 120 43 37];
            app.mveraLinkedin.ImageSource = 'linkedin_icon.svg';

            % Create chueteResearch
            app.chueteResearch = uiimage(app.AboutUIFigure);
            app.chueteResearch.ImageClickedFcn = createCallbackFcn(app, @chueteResearchClicked, true);
            app.chueteResearch.Position = [506 120 43 37];
            app.chueteResearch.ImageSource = 'researchgate_icon.svg';

            % Create chueteLinkedin
            app.chueteLinkedin = uiimage(app.AboutUIFigure);
            app.chueteLinkedin.ImageClickedFcn = createCallbackFcn(app, @chueteLinkedinClicked, true);
            app.chueteLinkedin.Position = [555 120 43 37];
            app.chueteLinkedin.ImageSource = 'linkedin_icon.svg';

            % Create acuadraImage
            app.acuadraImage = uiimage(app.AboutUIFigure);
            app.acuadraImage.ImageClickedFcn = createCallbackFcn(app, @acuadraWebClicked, true);
            app.acuadraImage.Position = [61 218 130 148];
            app.acuadraImage.ImageSource = 'acuadra_image.svg';

            % Create mveraImage
            app.mveraImage = uiimage(app.AboutUIFigure);
            app.mveraImage.ImageClickedFcn = createCallbackFcn(app, @mveraWebClicked, true);
            app.mveraImage.Position = [272 218 130 148];
            app.mveraImage.ImageSource = 'mvera_image.svg';

            % Create chueteImage
            app.chueteImage = uiimage(app.AboutUIFigure);
            app.chueteImage.ImageClickedFcn = createCallbackFcn(app, @chueteWebClicked, true);
            app.chueteImage.Position = [487 218 130 148];
            app.chueteImage.ImageSource = 'chuete_image.svg';

            % Create fluidosImage
            app.fluidosImage = uiimage(app.AboutUIFigure);
            app.fluidosImage.ImageClickedFcn = createCallbackFcn(app, @fluidosWebClicked, true);
            app.fluidosImage.Position = [16 26 464 51];
            app.fluidosImage.ImageSource = 'fluids_uc3m_icon.svg';

            % Show the figure after all components are created
            app.AboutUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = about

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.AboutUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.AboutUIFigure)
        end
    end
end
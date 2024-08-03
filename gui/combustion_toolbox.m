classdef combustion_toolbox < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        NewMenu                         matlab.ui.container.Menu
        OpenMenu                        matlab.ui.container.Menu
        SaveasMenu                      matlab.ui.container.Menu
        xlsMenu                         matlab.ui.container.Menu
        matMenu                         matlab.ui.container.Menu
        SnapshotMenu                    matlab.ui.container.Menu
        PreferencesMenu                 matlab.ui.container.Menu
        HelpMenu                        matlab.ui.container.Menu
        DocumentationMenu               matlab.ui.container.Menu
        TutorialMenu                    matlab.ui.container.Menu
        ExamplesMenu                    matlab.ui.container.Menu
        ValidationsMenu                 matlab.ui.container.Menu
        DiscussionsMenu                 matlab.ui.container.Menu
        CheckforupdatesMenu             matlab.ui.container.Menu
        LicenseMenu                     matlab.ui.container.Menu
        SendfeedbackMenu                matlab.ui.container.Menu
        AboutMenu                       matlab.ui.container.Menu
        Console_text                    matlab.ui.control.TextArea
        Lamp                            matlab.ui.control.Lamp
        Label_Console                   matlab.ui.control.Label
        Console                         matlab.ui.control.TextArea
        Tab_lateral_bar                 matlab.ui.container.TabGroup
        SetupTab                        matlab.ui.container.Tab
        TabGroup2                       matlab.ui.container.TabGroup
        InputsTab                       matlab.ui.container.Tab
        SelectProblemTypePanel          matlab.ui.container.Panel
        FLAG_IAC                        matlab.ui.control.CheckBox
        IonizedspeciesCheckBox          matlab.ui.control.CheckBox
        FrozenchemistryCheckBox         matlab.ui.control.CheckBox
        ProblemType                     matlab.ui.control.DropDown
        Calculate                       matlab.ui.control.Button
        Clear                           matlab.ui.control.Button
        AdditionalconstraintsPanel      matlab.ui.container.Panel
        PP5                             matlab.ui.control.EditField
        text_RP5                        matlab.ui.control.Label
        PR5                             matlab.ui.control.EditField
        text_P2                         matlab.ui.control.Label
        PP3                             matlab.ui.control.EditField
        PP4                             matlab.ui.control.EditField
        text_RP4                        matlab.ui.control.Label
        PR4                             matlab.ui.control.EditField
        text_RP3                        matlab.ui.control.Label
        PR3                             matlab.ui.control.EditField
        text_RP                         matlab.ui.control.Label
        text_R2                         matlab.ui.control.Label
        DefinestateofreactantsandproductsPanel  matlab.ui.container.Panel
        text_RP2_2                      matlab.ui.control.Label
        text_RP1_2                      matlab.ui.control.Label
        text_RP2                        matlab.ui.control.Label
        PR2                             matlab.ui.control.EditField
        text_P1                         matlab.ui.control.Label
        text_R1                         matlab.ui.control.Label
        text_RP1                        matlab.ui.control.Label
        PR1                             matlab.ui.control.EditField
        PP1                             matlab.ui.control.EditField
        PP2                             matlab.ui.control.EditField
        DefinereactantsandspeciestobeconsideredPanel  matlab.ui.container.Panel
        PeriodicTable_R                 matlab.ui.control.Button
        edit_F                          matlab.ui.control.NumericEditField
        ProductsPanel                   matlab.ui.container.Panel
        PeriodicTable_P                 matlab.ui.control.Button
        Products                        matlab.ui.control.DropDown
        FuelEditFieldLabel              matlab.ui.control.Label
        text_phi                        matlab.ui.control.Label
        edit_phi                        matlab.ui.control.EditField
        edit_OF                         matlab.ui.control.NumericEditField
        OFEditFieldLabel                matlab.ui.control.Label
        ListofSpeciesPanel              matlab.ui.container.Panel
        listbox_Products                matlab.ui.control.ListBox
        ReactantsPanel                  matlab.ui.container.Panel
        Reactants                       matlab.ui.control.DropDown
        UITable_R                       matlab.ui.control.Table
        QuicksettingsTab                matlab.ui.container.Tab
        ReportPanel                     matlab.ui.container.Panel
        PrintresultsCheckBox            matlab.ui.control.CheckBox
        Report_type                     matlab.ui.control.DropDown
        TypeDropDownLabel               matlab.ui.control.Label
        OthersPanel                     matlab.ui.container.Panel
        IdealAirCheckBox                matlab.ui.control.CheckBox
        TuningparametersPanel           matlab.ui.container.Panel
        MaxiterationsSDEditField        matlab.ui.control.NumericEditField
        MaxiterationsSDEditFieldLabel   matlab.ui.control.Label
        RFMT0EditField                  matlab.ui.control.NumericEditField
        RFMT0EditFieldLabel             matlab.ui.control.Label
        RFMT0_LEditField                matlab.ui.control.NumericEditField
        RFMT0_LEditFieldLabel           matlab.ui.control.Label
        RFMT0_REditField                matlab.ui.control.NumericEditField
        RFMT0_REditFieldLabel           matlab.ui.control.Label
        MaxiterationsRFMEditField       matlab.ui.control.NumericEditField
        MaxiterationsRFMEditFieldLabel  matlab.ui.control.Label
        TolerancesPanel                 matlab.ui.container.Panel
        ShocksandDetonationsEditField   matlab.ui.control.NumericEditField
        ShocksandDetonationsEditFieldLabel  matlab.ui.control.Label
        RootFindingMethodEditField      matlab.ui.control.NumericEditField
        RootFindingMethodEditFieldLabel  matlab.ui.control.Label
        DisplaySpeciesEditField         matlab.ui.control.NumericEditField
        DisplaySpeciesEditFieldLabel    matlab.ui.control.Label
        TraceoptionEditField            matlab.ui.control.NumericEditField
        TraceoptionEditFieldLabel       matlab.ui.control.Label
        TabGroup4                       matlab.ui.container.TabGroup
        SelectionofspeciesTab           matlab.ui.container.Tab
        PeriodicTable_P_2               matlab.ui.control.Button
        listbox_LS                      matlab.ui.control.ListBox
        text_LS                         matlab.ui.control.Label
        RemoveButton1                   matlab.ui.control.Button
        AddButton1                      matlab.ui.control.Button
        listbox_LS_DB                   matlab.ui.control.ListBox
        text_LS_DB                      matlab.ui.control.Label
        edit_seeker_1                   matlab.ui.control.EditField
        text_list_species               matlab.ui.control.Label
        DisplayspeciesTab               matlab.ui.container.Tab
        text_LS_display                 matlab.ui.control.Label
        listbox_LS_display              matlab.ui.control.ListBox
        RemoveButton2                   matlab.ui.control.Button
        AddButton2                      matlab.ui.control.Button
        listbox_LS_2                    matlab.ui.control.ListBox
        text_LS_2                       matlab.ui.control.Label
        edit_seeker_2                   matlab.ui.control.EditField
        text_list_species_2             matlab.ui.control.Label
        ResultsTab                      matlab.ui.container.Tab
        Tree                            matlab.ui.container.Tree
        Node_Results                    matlab.ui.container.TreeNode
        Tab_results                     matlab.ui.container.TabGroup
        ParametersTab                   matlab.ui.container.Tab
        Panel_parameters                matlab.ui.container.Panel
        text_theta_2                    matlab.ui.control.NumericEditField
        text_beta_2                     matlab.ui.control.NumericEditField
        text_beta_min_2                 matlab.ui.control.NumericEditField
        text_theta                      matlab.ui.control.Label
        text_beta                       matlab.ui.control.Label
        text_beta_min                   matlab.ui.control.Label
        text_Isp_2                      matlab.ui.control.NumericEditField
        text_Ivac_2                     matlab.ui.control.NumericEditField
        text_Cstar_2                    matlab.ui.control.NumericEditField
        EpsilonmolesLabel_2             matlab.ui.control.Label
        text_error_problem              matlab.ui.control.NumericEditField
        Panel_extra_5                   matlab.ui.container.Panel
        text_Isp_5                      matlab.ui.control.NumericEditField
        text_Ivac_5                     matlab.ui.control.NumericEditField
        text_Cstar_5                    matlab.ui.control.NumericEditField
        text_Aratio_5                   matlab.ui.control.NumericEditField
        text_T_5                        matlab.ui.control.NumericEditField
        text_Products_5                 matlab.ui.control.Label
        text_p_5                        matlab.ui.control.NumericEditField
        text_r_5                        matlab.ui.control.NumericEditField
        text_h_5                        matlab.ui.control.NumericEditField
        text_cp_5                       matlab.ui.control.NumericEditField
        text_W_5                        matlab.ui.control.NumericEditField
        text_sound_5                    matlab.ui.control.NumericEditField
        text_e_5                        matlab.ui.control.NumericEditField
        text_u_5                        matlab.ui.control.NumericEditField
        text_M_5                        matlab.ui.control.NumericEditField
        text_gamma_5                    matlab.ui.control.NumericEditField
        text_s_5                        matlab.ui.control.NumericEditField
        Panel_extra_4                   matlab.ui.container.Panel
        text_Isp_4                      matlab.ui.control.NumericEditField
        text_Ivac_4                     matlab.ui.control.NumericEditField
        text_Cstar_4                    matlab.ui.control.NumericEditField
        text_Aratio_4                   matlab.ui.control.NumericEditField
        text_T_4                        matlab.ui.control.NumericEditField
        text_Products_4                 matlab.ui.control.Label
        text_p_4                        matlab.ui.control.NumericEditField
        text_r_4                        matlab.ui.control.NumericEditField
        text_h_4                        matlab.ui.control.NumericEditField
        text_cp_4                       matlab.ui.control.NumericEditField
        text_W_4                        matlab.ui.control.NumericEditField
        text_sound_4                    matlab.ui.control.NumericEditField
        text_e_4                        matlab.ui.control.NumericEditField
        text_u_4                        matlab.ui.control.NumericEditField
        text_M_4                        matlab.ui.control.NumericEditField
        text_gamma_4                    matlab.ui.control.NumericEditField
        text_s_4                        matlab.ui.control.NumericEditField
        Panel_extra_3                   matlab.ui.container.Panel
        text_Isp_3                      matlab.ui.control.NumericEditField
        text_Ivac_3                     matlab.ui.control.NumericEditField
        text_Cstar_3                    matlab.ui.control.NumericEditField
        text_Aratio_3                   matlab.ui.control.NumericEditField
        text_T_3                        matlab.ui.control.NumericEditField
        text_Products_3                 matlab.ui.control.Label
        text_p_3                        matlab.ui.control.NumericEditField
        text_r_3                        matlab.ui.control.NumericEditField
        text_h_3                        matlab.ui.control.NumericEditField
        text_cp_3                       matlab.ui.control.NumericEditField
        text_W_3                        matlab.ui.control.NumericEditField
        text_sound_3                    matlab.ui.control.NumericEditField
        text_e_3                        matlab.ui.control.NumericEditField
        text_u_3                        matlab.ui.control.NumericEditField
        text_M_3                        matlab.ui.control.NumericEditField
        text_gamma_3                    matlab.ui.control.NumericEditField
        text_s_3                        matlab.ui.control.NumericEditField
        text_Isp                        matlab.ui.control.Label
        text_Ivac                       matlab.ui.control.Label
        text_Cstar                      matlab.ui.control.Label
        text_Aratio_2                   matlab.ui.control.NumericEditField
        text_Aratio                     matlab.ui.control.Label
        text_M_1                        matlab.ui.control.NumericEditField
        text_u_1                        matlab.ui.control.NumericEditField
        text_T_2                        matlab.ui.control.NumericEditField
        text_p_1                        matlab.ui.control.NumericEditField
        text_T_1                        matlab.ui.control.NumericEditField
        text_T                          matlab.ui.control.Label
        text_Reactans                   matlab.ui.control.Label
        text_Products                   matlab.ui.control.Label
        text_p_2                        matlab.ui.control.NumericEditField
        text_p                          matlab.ui.control.Label
        text_r_2                        matlab.ui.control.NumericEditField
        text_r_1                        matlab.ui.control.NumericEditField
        text_r                          matlab.ui.control.Label
        text_h_2                        matlab.ui.control.NumericEditField
        text_h_1                        matlab.ui.control.NumericEditField
        text_h                          matlab.ui.control.Label
        text_cp_2                       matlab.ui.control.NumericEditField
        text_cp_1                       matlab.ui.control.NumericEditField
        text_cp                         matlab.ui.control.Label
        text_W_2                        matlab.ui.control.NumericEditField
        text_W_1                        matlab.ui.control.NumericEditField
        text_W                          matlab.ui.control.Label
        text_sound_2                    matlab.ui.control.NumericEditField
        text_sound_1                    matlab.ui.control.NumericEditField
        text_sound                      matlab.ui.control.Label
        text_e_2                        matlab.ui.control.NumericEditField
        text_e_1                        matlab.ui.control.NumericEditField
        text_e                          matlab.ui.control.Label
        text_u_2                        matlab.ui.control.NumericEditField
        text_u                          matlab.ui.control.Label
        text_M_2                        matlab.ui.control.NumericEditField
        text_M                          matlab.ui.control.Label
        text_gamma_2                    matlab.ui.control.NumericEditField
        text_gamma_1                    matlab.ui.control.NumericEditField
        text_gamma                      matlab.ui.control.Label
        text_s_2                        matlab.ui.control.NumericEditField
        text_s_1                        matlab.ui.control.NumericEditField
        text_s                          matlab.ui.control.Label
        text_phi_3                      matlab.ui.control.Label
        edit_phi3                       matlab.ui.control.EditField
        MixturecompositionTab           matlab.ui.container.Tab
        text_phi_2                      matlab.ui.control.Label
        edit_phi2                       matlab.ui.control.EditField
        text_R3                         matlab.ui.control.Label
        text_P3                         matlab.ui.control.Label
        UITable_P                       matlab.ui.control.Table
        UITable_R2                      matlab.ui.control.Table
        CustomFiguresTab                matlab.ui.container.Tab
        DefaultsettingsCheckBox         matlab.ui.control.CheckBox
        figure_settings                 matlab.ui.control.Button
        figure_size                     matlab.ui.control.Button
        figure_clear                    matlab.ui.control.Button
        figure_plot                     matlab.ui.control.Button
        Tree_variable_y                 matlab.ui.container.CheckBoxTree
        Variable_y                      matlab.ui.container.TreeNode
        Tree_variable_x                 matlab.ui.container.CheckBoxTree
        Variable_x                      matlab.ui.container.TreeNode
        Tree_mixtures                   matlab.ui.container.CheckBoxTree
        Mixtures                        matlab.ui.container.TreeNode
        UIAxes                          matlab.ui.control.UIAxes
        ContextMenu_CommandWindow       matlab.ui.container.ContextMenu
        MaximizeMenu                    matlab.ui.container.Menu
        MinimizeMenu                    matlab.ui.container.Menu
        ClearCommandWindowMenu          matlab.ui.container.Menu
        ContextMenu_UITree              matlab.ui.container.ContextMenu
        RemoveMenu                      matlab.ui.container.Menu
        ContextMenu_UIAxes              matlab.ui.container.ContextMenu
        property_inspector_menu         matlab.ui.container.Menu
        Menu                            matlab.ui.container.Menu
    end

    properties (Access = public)
        constants         % Constants object
        database          % Class inhereted from Database superclass
        chemicalSystem    % Chemical system object
        mixture           % Mixture object
        equilibriumSolver % EquilibriumSolver object
        shockSolver       % ShockSolver object
        detonationSolver  % DetonationSolver object
        rocketSolver      % RocketSolver object
        plotConfig        % PlotConfig object
        export            % Export object
        displaySpecies    % Checmial species to be shown in plotComposition
        fig               % Auxiliary figure
        default           % Struct with default values of some components in the GUI
        N_flags           % Number of flags active
        PR1_vector        % Condition Reactants 1
        PR2_vector        % Condition Reactants 2
        PR3_vector        % Condition Reactants 3
        PP1_vector        % Condition Products 1
        PP2_vector        % Condition Products 2
        PR1_var_name      % Variable name for PR1
        PR2_var_name      % Variable name for PR2
        PR3_var_name      % Variable name for PR3
        PP1_var_name      % Variable name for PP1
        PP2_var_name      % Variable name for PP2
        flag_PR1          % FLAG for PR1: true-> vector
        flag_PR2          % FLAG for PR2: true-> vector
        flag_PR3          % FLAG for PR3: true-> vector
        flag_PP1          % FLAG for PP1: true-> vector
        flag_PP2          % FLAG for PP2: true-> vector
        flag_phi          % FLAG for phi: true-> vector
        indexFuel         % Index position Fuel species
        indexOxidizer     % Index position Oxidizer species
        indexInert        % Index position Inert species
        LS                % List of species considered (reactants + products)
        LS_products       % List of species considered as products
        color_splash       = [0.5098, 0.6039, 0.6745]; % Font color splash
        color_lamp_nothing = [0.8000, 0.8000, 0.8000]; % Lamp color (rgb): nothing to report
        color_lamp_working = [0.9961, 0.9804, 0.8314]; % Lamp color (rgb): working
        color_lamp_done    = [0.5608, 0.7255, 0.6588]; % Lamp color (rgb): done
        color_lamp_error   = [0.6400, 0.0800, 0.1800]; % Lamp color (rgb): error
        temp_results       % Temporal variable that contains the last parametric study
        current_history    % Current history of commands
        temp_index         % Temporal index to get current position of command history
        dynamic_components % Struct with all the dynamic components
        welcome_message = 'Welcome to Combustion Toolbox %s --- A MATLAB-GUI based open-source tool for solving gaseous combustion problems.';
        maxRelativeError = 2e-2; % Relative error threshold to change color (2 %)
    end
    
    properties (Dependent)
        NS_products % Number of product species (computations)
        NS_display  % Number of display species (plots)
        LS_reactants % List of reactants
    end

    methods (Access = public)
        % Value changed function: Reactants
        function public_ReactantsValueChanged(app, event)
            gui_ReactantsValueChanged(app, event);
            % Compute pre-shock velocity (only for shocks)
            gui_compute_mach_or_velocity(app, 'Mach');
            % Update Listbox and Display species (extended settings)
            public_ProductsValueChanged(app);
        end

        % Value changed function: Products
        function public_ProductsValueChanged(app)
            % Update Listbox (extended settings)
            app.listbox_LS.Items = app.listbox_Products.Items;
            % * Update Title with the number of items contained in the box
            app.ListofSpeciesPanel.Title = sprintf('List of Species - %d', app.NS_products);
            % Update Display species (extended settings)
            app.listbox_LS_2.Items = app.listbox_LS.Items;
            % Set Display species to all (extended settings)
            app.listbox_LS_display.Items = app.listbox_LS_2.Items;
            app.displaySpecies = app.listbox_LS_display.Items;
            % Update Text with the number of items contained in the box
            app.text_LS.Text = sprintf('List of Species - %d', app.NS_products);
            % Update Text with the number of items contained in the box
            app.text_LS_2.Text = sprintf('List of Species - %d', app.NS_products);
            % Update Text with the number of items contained in the box
            app.text_LS_display.Text = sprintf('Display Species - %d', app.NS_display);
        end

        function public_FLAG_IACValueChanged(app)
            if ~app.FLAG_IAC.Value
                app.text_P1.Text = 'FAC';
                app.text_P1.Visible = 'on';
                app.text_RP1_2.Visible = 'on'; app.text_RP2_2.Visible = 'off';
                app.PP1.Visible = 'on'; app.PP2.Visible = 'off';
                app.PP1.Value = ''; app.PP2.Value = '';
                app.Panel_extra_5.Visible = 'on';
                app.text_Products.Text = 'Injector';
                app.text_Products_3.Text = 'Outlet Chamber';
                app.text_Products_4.Text = 'Throat';
                app.text_Products_5.Text = 'Exit';
            else
                app.text_P1.Text = 'Products';
                app.text_P1.Visible = 'off';
                app.text_RP1_2.Visible = 'off'; app.text_RP2_2.Visible = 'off';
                app.PP1.Visible = 'off'; app.PP2.Visible = 'off';
                app.PP1.Value = '2500'; app.PP2.Value = '1';
                app.Panel_extra_5.Visible = 'off';
                app.text_Products.Text = 'Outlet Chamber';
                app.text_Products_3.Text = 'Throat';
                app.text_Products_4.Text = 'Exit';
                app.text_Products_5.Text = 'Exit';
            end
        end

        function [current_history, temp_index] = public_get_current_history(app)
            % Get current history
            current_history = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
            % Get length history
            temp_index = length(current_history);
            % Assign values
            app.current_history = current_history;
            app.temp_index = temp_index;
        end
        
        function public_ClearButtonPushed(app, event)
            ClearButtonPushed(app, event)
        end
        
        function public_xlsMenuSelected(app, event)
            xlsMenuSelected(app, event);
        end

        function public_matMenuSelected(app, event)
            matMenuSelected(app, event);
        end

    end

    methods
        function value = get.NS_products(app)
            value = length(app.listbox_Products.Items);
        end

        function value = get.NS_display(app)
            value = length(app.displaySpecies);
        end

        function value = get.LS_reactants(app)
            try
                value = app.UITable_R.Data(:,1)';
            catch
                value = [];
            end
            
        end

    end

    methods (Access = private)
        
        function create_components(app)
            % Create additional components
            
            % Create toolbar background
            app.dynamic_components.toolbar_background = uiimage(app.UIFigure, 'Position', [485 746 180 20], 'ImageSource', 'toolbar_background.svg');
            % Add icons
            app.dynamic_components.toolbar_button1 = uiimage(app.UIFigure, 'Position', [520 746 20 20], 'ImageSource', 'icon_new.svg', 'Tooltip', 'New (Ctrl + N)', 'ImageClickedFcn', createCallbackFcn(app, @ClearButtonPushed, true));
            app.dynamic_components.toolbar_button2 = uiimage(app.UIFigure, 'Position', [542 746 20 20], 'ImageSource', 'icon_save.svg', 'Tooltip', 'Save (Ctrl + S)', 'ImageClickedFcn', createCallbackFcn(app, @xlsMenuSelected, true));
            app.dynamic_components.toolbar_button3 = uiimage(app.UIFigure, 'Position', [564 746 20 20], 'ImageSource', 'icon_play.svg', 'Tooltip', 'Calculate (F5)', 'ImageClickedFcn', createCallbackFcn(app, @CalculateButtonPushed, true));
            app.dynamic_components.toolbar_button4 = uiimage(app.UIFigure, 'Position', [586 746 20 20], 'ImageSource', 'icon_help.svg', 'Tooltip', 'Help (F1)', 'ImageClickedFcn', @(~,~) system('start https://combustion-toolbox-website.readthedocs.io/_/downloads/en/latest/pdf/'));
            app.dynamic_components.toolbar_button5 = uiimage(app.UIFigure, 'Position', [608 746 20 20], 'ImageSource', 'icon_github_CT.svg', 'Tooltip', 'GitHub', 'ImageClickedFcn', @(~,~) system('start https://github.com/AlbertoCuadra/combustion_toolbox'));
        end

    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Load figure in the background
            app.UIFigure.Visible = 'off';
            
            % Splash
            try
                FLAG_SPLASH = true;
                splash_obj = gui_display_splash('app', app);
            catch
                FLAG_SPLASH = false;
            end

            % Add additional components
            % create_components(app); % Next release
            % Get screen position
            % position = splash_obj.ScreenPosition;
            position = combustiontoolbox.utils.display.getMonitorPositionsMATLAB();
            % Centering app
            app.UIFigure.Position(1) = position(1) + (position(3) - app.UIFigure.Position(3))/2;
            app.UIFigure.Position(2) = position(2) + (position(4) - app.UIFigure.Position(4))/2;
            % Check screen size
            if position(4) <= app.UIFigure.Position(4)
                app.UIFigure.Position(4) = 0.9 * position(4);
                app.UIFigure.Resize = true;
                app.UIFigure.Scrollable = true;
            end
            % Get Constants
            app.constants = combustiontoolbox.common.Constants;
            % Get Nasa's database
            app.database = combustiontoolbox.databases.NasaDatabase();
            % Initialize chemical system
            app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database);
            % Initialize mixture
            app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);
            % Initialize plotConfig object
            app.plotConfig = combustiontoolbox.utils.display.PlotConfig;
            % Initialize solvers
            app.equilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('plotConfig', app.plotConfig);
            app.shockSolver = combustiontoolbox.shockdetonation.ShockSolver('plotConfig', app.plotConfig, 'equilibriumSolver', app.equilibriumSolver);
            app.detonationSolver = combustiontoolbox.shockdetonation.DetonationSolver('plotConfig', app.plotConfig, 'equilibriumSolver', app.equilibriumSolver);
            app.rocketSolver = combustiontoolbox.rocket.RocketSolver('plotConfig', app.plotConfig, 'equilibriumSolver', app.equilibriumSolver);
            % Initialize export object
            app.export = combustiontoolbox.utils.Export();
            % Initialize List box species DataBase master
            app.listbox_LS_DB.Items = app.database.listSpecies;
            % Initialize table data
            app.UITable_R.ColumnFormat = {[] [] [] {'Fuel', 'Oxidizer', 'Inert'} []};
            % Initialize for TP
            app.PR1_var_name = 'TR'; app.PR2_var_name = 'pR';
            app.PP1_var_name = 'TP'; app.PP2_var_name = 'pP';
            % Get tolerances
            % app = gui_get_tolerances(app); % Remove
            % Save default state GUI
            app.default.Panel_parameters.Position = app.Panel_parameters.Position;
            app.default.data_UITable_R = app.UITable_R.Data;
%             app.UITable_R.ColumnFormat(2:3)  = {'shortE'};
%             app.UITable_R2.ColumnFormat(2:3) = {'shortE'};
            app.UITable_P.ColumnFormat(2:3)  = {'shortE'};
            % Update welcome message in the command window
            app.Console_text.Value = sprintf(app.welcome_message, app.constants.release);
            % Get current history
            public_get_current_history(app);
            % Make app visible
            app.UIFigure.Visible = 'on';
            % Delete splash
            if FLAG_SPLASH
                delete(splash_obj);
            end
        end

        % Menu selected function: AboutMenu
        function AboutMenuSelected(app, event)
            uiabout();
        end

        % Value changed function: Console
        function ConsoleValueChanged(app, event)
            gui_ConsoleValueChanged(app, event);
        end

        % Menu selected function: SnapshotMenu
        function SnapshotMenuSelected(app, event)
            gui_SnapshotMenuSelected(app.UIFigure);
        end

        % Menu selected function: MaximizeMenu
        function MaximizeMenuSelected(app, event)
            app.Console_text.Position = [78, 30, 554, 462];
        end

        % Menu selected function: MinimizeMenu
        function MinimizeMenuSelected(app, event)
            app.Console_text.Position = [78, 30, 554 57];
        end

        % Menu selected function: ClearCommandWindowMenu
        function ClearCommandWindowMenuSelected(app, event)
            app.Console_text.Value = '';
        end

        % Button pushed function: Calculate
        function CalculateButtonPushed(app, event)
%             profile on
            gui_CalculateButtonPushed(app, event);
%             profile viewer
        end

        % Menu selected function: SendfeedbackMenu
        function SendfeedbackMenuSelected(app, event)
            run('uifeedback');
        end

        % Cell edit callback: UITable_R
        function UITable_RCellEdit(app, event)
            gui_UITable_RCellEdit(app, event);
        end

        % Value changed function: edit_phi
        function edit_phiValueChanged(app, event)
            gui_edit_phiValueChanged(app, event);
        end

        % Value changed function: PR1
        function PR1ValueChanged(app, event)
            % Update equilibrium temperature of reactants (GUI)            
            if ~isempty(app.UITable_R.Data)
                listSpecies = app.UITable_R.Data(:, 1);
                temperature = {gui_get_prop(app.PR1.Value)};
                temperature = {temperature{1}(1)};
                numSpecies = length(app.UITable_R.Data(:, 1));
                [app, temperature, FLAG] = gui_check_temperature_reactants(app, listSpecies, temperature, numSpecies);
                
                if FLAG
                    app.UITable_R.Data(:, 5) = temperature;
                    message = sprintf('The species selected as reactants can only be evaluated at its defined temperature.\nThe temperature shown as the temperature of the reactant is a ficticious value! The species are unmixed.');
                    uialert(app.UIFigure, message, 'Warning', 'Icon', 'warning');
                else
                    app.UITable_R.Data(:, 5) = repmat(temperature(1), [1, numSpecies]);
                end

            end

            % Compute pre-shock velocity (only for shocks)
            gui_compute_mach_or_velocity(app, 'Mach');
        end

        % Menu selected function: DiscussionsMenu
        function DiscussionsMenuSelected(app, event)
            % Open default web browser and redirect to Combustion Toolbox
            % discussion forum

            % Import packages
            import combustiontoolbox.utils.SystemUtils
        
            % Definitions
            url = SystemUtils.url.discussions;
            
            % Open website
            SystemUtils.openWebsite(url)
        end

        % Value changed function: ProblemType
        function ProblemTypeValueChanged(app, event)
            % Clear GUI results tab (except UITree) and update GUI items for the problem selected
            gui_ProblemTypeValueChanged(app);
        end

        % Value changed function: PP2
        function PP2ValueChanged(app, event)
            if app.ProblemType.Value(2) == 'P' &&  app.ProblemType.Value(1) ~= 'S'
                app.PR2.Value = event.Value;
            end
        end

        % Value changed function: PR2
        function PR2ValueChanged(app, event)
            if app.ProblemType.Value(2) == 'P' &&  app.ProblemType.Value(1) ~= 'S'
                app.PP2.Value = event.Value;
            end
            % Compute pre-shock velocity (only for shocks)
            gui_compute_mach_or_velocity(app, 'Mach');
        end

        % Value changed function: FrozenchemistryCheckBox
        function FrozenchemistryCheckBoxValueChanged(app, event)
            % Set frozen chemistry
            gui_FrozenchemistryCheckBoxValueChanged(app)
        end

        % Menu selected function: DocumentationMenu
        function DocumentationMenuSelected(app, event)
            combustiontoolbox.utils.SystemUtils.websiteCT;
        end

        % Button pushed function: Clear
        function ClearButtonPushed(app, event)
            % Clear command window
            ClearCommandWindowMenuSelected(app, event);
            % Clear axes
            cla(app.UIAxes);
            % Delete nodes uitree
            delete(app.Node_Results.Children);
            delete(app.Mixtures.Children);
            delete(app.Variable_x.Children);
            delete(app.Variable_y.Children);
            % Set color lamp to nothing
            app.Lamp.Color = app.color_lamp_nothing;
        end

        % Selection changed function: Tree
        function TreeSelectionChanged(app, event)
            selectedNodes = app.Tree.SelectedNodes;
            gui_update_from_uitree(app, selectedNodes(1));
        end

        % Value changed function: PR4
        function PR4ValueChanged(app, event)
            % Compute pre-shock velocity (only for shocks)
            gui_compute_mach_or_velocity(app, 'Mach');
        end

        % Value changed function: PR3
        function PR3ValueChanged(app, event)
            % Compute pre-shock Mach number (only for shocks)
            gui_compute_mach_or_velocity(app, 'velocity');
        end

        % Menu selected function: xlsMenu
        function xlsMenuSelected(app, event)
            % Extract last results
            mixArray1 = [app.temp_results.mix1];
            mixArray2 = [app.temp_results.mix2];
            
            % Store default format
            temp = app.export.format;

            % Export results
            app.export.format = '.xls';
            app.export.export(mixArray1, mixArray2);

            % Recover default format
            app.export.format = temp;
        end

        % Menu selected function: matMenu
        function matMenuSelected(app, event)
            % Extract last results
            mixArray1 = [app.temp_results.mix1];
            mixArray2 = [app.temp_results.mix2];
            
            % Store default format
            temp = app.export.format;

            % Export results
            app.export.format = '.mat';
            app.export.export(mixArray1, mixArray2);

            % Recover default format
            app.export.format = temp;
        end

        % Button pushed function: AddButton1
        function AddButton1Pushed(app, event)
            app.listbox_LS.Items = gui_value2list(app, app.listbox_LS_DB.Value, app.listbox_LS.Items, 'add');
            % Update Listbox (inputs)
            app.listbox_Products.Items = app.listbox_LS.Items;
            % Update counters
            public_ProductsValueChanged(app);
        end

        % Button pushed function: RemoveButton1
        function RemoveButton1Pushed(app, event)
            app.listbox_LS.Items = gui_value2list(app, app.listbox_LS.Value, app.listbox_LS.Items, 'remove');
            % Update Listbox (inputs)
            app.listbox_Products.Items = app.listbox_LS.Items;
            % Update counters
            public_ProductsValueChanged(app);
        end

        % Menu selected function: ValidationsMenu
        function ValidationsMenuSelected(app, event)
            uivalidations();
        end

        % Menu selected function: RemoveMenu
        function RemoveMenuSelected(app, event)
            delete(app.Tree.SelectedNodes);
        end

        % Menu selected function: LicenseMenu
        function LicenseMenuSelected(app, event)
            uilicense();
        end

        % Button pushed function: figure_clear
        function figure_clearButtonPushed(app, event)
            ClearButtonPushed(app, event);
        end

        % Button pushed function: figure_plot
        function figure_plotButtonPushed(app, event)
            gui_plot_custom_figures(app);
        end

        % Button pushed function: figure_size
        function figure_sizeButtonPushed(app, event)
            if contains(app.figure_size.Text, 'Max')
                % Save default values
                app.default.UIFigure_position = app.UIFigure.Position;
                app.default.UIAxes_position = app.UIAxes.Position;
                app.default.Tab_results_position = app.Tab_results.Position;
                % Maximize UIFigure (width)
                app.UIFigure.Position([1,3]) = [0.5 * app.UIFigure.Position(1), 2 * app.UIFigure.Position(3)];
                % Delete UIAxes
                delete(app.UIAxes);
                % Create new UIAxes
                app.UIAxes = uiaxes(app.UIFigure, 'Position', [655,9,609,719]);
                % Change name button
                app.figure_size.Text = 'Minimize';
                % Replot
                figure_plotButtonPushed(app, event);
            else
                app.UIFigure.Position = app.default.UIFigure_position;
                % Delete UIAxes
                delete(app.UIAxes);
                % Create new UIAxes
                app.UIAxes = uiaxes(app.CustomFiguresTab, 'Position', app.default.UIAxes_position);
                % Change name button
                app.figure_size.Text = 'Maximize';
                % Replot
                figure_plotButtonPushed(app, event);
            end
        end

        % Menu selected function: CheckforupdatesMenu
        function CheckforupdatesMenuSelected(app, event)
            [~, message] = combustiontoolbox.utils.checkUpdate(app.UIFigure);
            app.Console_text.Value = message;
        end

        % Button pushed function: PeriodicTable_R
        function PeriodicTable_RButtonPushed(app, event)
            uielements(app, event);
        end

        % Button pushed function: PeriodicTable_P, PeriodicTable_P_2
        function PeriodicTable_PButtonPushed(app, event)
            uielements(app, event);
        end

        % Value changed function: Reactants
        function ReactantsValueChanged(app, event)
            public_ReactantsValueChanged(app, event);
        end

        % Value changed function: Products
        function ProductsValueChanged(app, event)
            % Update List of species considered as Products
            gui_ProductsValueChanged(app);
            % Update Listbox (extended settings)
            public_ProductsValueChanged(app);
        end

        % Value changing function: edit_seeker_1
        function edit_seeker_1ValueChanging(app, event)
            seek_value = gui_seeker_value(app, event, app.database.listSpecies);
            % Update Listbox (inputs)
            if isempty(event.Value)
                app.listbox_LS_DB.Items = app.database.listSpecies;
                return
            end
            
            app.listbox_LS_DB.Items = seek_value;
        end

        % Value changing function: edit_seeker_2
        function edit_seeker_2ValueChanging(app, event)
            seek_value1 = gui_seeker_value(app, event, app.listbox_LS.Items);
            seek_value2 = gui_seeker_value(app, event, app.displaySpecies);
            % Update Listbox (inputs)
            if isempty(event.Value)
                app.listbox_LS_2.Items = app.listbox_LS.Items;
                app.listbox_LS_display.Items = app.displaySpecies;
                return
            end
            
            app.listbox_LS_2.Items = seek_value1;
            app.listbox_LS_display.Items = seek_value2;
        end

        % Button pushed function: AddButton2
        function AddButton2Pushed(app, event)
            % Add species from LS to display species
            app.listbox_LS_display.Items = gui_value2list(app, app.listbox_LS_2.Value, app.listbox_LS_display.Items, 'add');
            % Update displaySpecies
            app.displaySpecies = unique([app.displaySpecies, app.listbox_LS_display.Items]);
            % Update Text with the number of items contained in the box
            app.text_LS_display.Text = sprintf('Display Species - %d', app.NS_display);
        end

        % Button pushed function: RemoveButton2
        function RemoveButton2Pushed(app, event)
            % Remove species from display species
            app.displaySpecies = gui_value2list(app, app.listbox_LS_display.Value, app.displaySpecies, 'remove');
            % Update items listbox_LS_display
            app.listbox_LS_display.Items = setdiff(app.listbox_LS_display.Items, app.listbox_LS_display.Value);
            % Update Text with the number of items contained in the box
            app.text_LS_display.Text = sprintf('Display Species - %d', app.NS_display);
        end

        % Menu selected function: property_inspector_menu
        function property_inspector_menuSelected(app, event)
            % Set default settings to 'off'
            app.DefaultsettingsCheckBox.Value = false;
            % Show property inspector
            inspect(app.UIAxes)
        end

        % Button pushed function: figure_settings
        function figure_settingsButtonPushed(app, event)
            % Set default settings to 'off'
            app.DefaultsettingsCheckBox.Value = false;
            % Show property inspector
            inspect(app.UIAxes)
        end

        % Value changed function: FLAG_IAC
        function FLAG_IACValueChanged(app, event)
            public_FLAG_IACValueChanged(app);            
        end

        % Value changing function: PR3
        function PR3ValueChanging(app, event)
            % Check only PR3 or PP3 not both (only for ROCKET)
            gui_keep_last_entry(app, app.PP3)
        end

        % Value changing function: PP3
        function PP3ValueChanging(app, event)
            % Check only PR4 or PP4 not both (only for ROCKET)
            gui_keep_last_entry(app, app.PR3)
        end

        % Callback function
        function PeriodicTable_displayButtonPushed(app, event)
            uielements(app, event);
        end

        % Menu selected function: PreferencesMenu
        function PreferencesMenuSelected(app, event)
            uipreferences(app);
        end

        % Value changed function: TraceoptionEditField
        function TraceoptionEditFieldValueChanged(app, event)
            app.equilibriumSolver.tolGibbs = app.TraceoptionEditField.Value;
        end

        % Value changed function: RootFindingMethodEditField
        function RootFindingMethodEditFieldValueChanged(app, event)
            app.equilibriumSolver.tol0 = app.RootFindingMethodEditField.Value;
        end

        % Value changed function: ShocksandDetonationsEditField
        function ShocksandDetonationsEditFieldValueChanged(app, event)
            app.shockSolver.tol0 = app.ShocksandDetonationsEditField.Value;
            app.detonationSolver.tol0 = app.ShocksandDetonationsEditField.Value;
        end

        % Value changed function: DisplaySpeciesEditField
        function DisplaySpeciesEditFieldValueChanged(app, event)
            app.plotConfig.mintolDisplay = app.DisplaySpeciesEditField.Value;
        end

        % Value changed function: MaxiterationsRFMEditField
        function MaxiterationsRFMEditFieldValueChanged(app, event)
            app.equilibriumSolver.itMax = app.MaxiterationsRFMEditField.Value;
        end

        % Value changed function: MaxiterationsSDEditField
        function MaxiterationsSDEditFieldValueChanged(app, event)
            app.shockSolver.itMax = app.MaxiterationsSDEditField.Value;
            app.detonationSolver.itMax = app.MaxiterationsSDEditField.Value;
        end

        % Value changed function: RFMT0_LEditField
        function RFMT0_LEditFieldValueChanged(app, event)
            app.equilibriumSolver.root_T0_l = app.RFMT0_LEditField.Value;
        end

        % Value changed function: RFMT0_REditField
        function RFMT0_REditFieldValueChanged(app, event)
            app.equilibriumSolver.root_T0_r = app.RFMT0_REditField.Value;
        end

        % Value changed function: RFMT0EditField
        function RFMT0EditFieldValueChanged(app, event)
            app.equilibriumSolver.root_T0 = app.RFMT0EditField.Value;
        end

        % Key press function: UIFigure
        function UIFigureKeyPress(app, event)
            % Definitions
            FLAG_ENTER = false;
            % Get key
            key = event.Key;
            try
                modifier = event.Modifier{:};
            catch
                modifier = '';
            end
            % Read action

            % Modifier used (shift + key press)
            if isequal(modifier,'shift')
                switch lower(key)
                    case 'f1'
                        event.Source.Tag = 'R';
                        PeriodicTable_RButtonPushed(app, event);
                    case 'f2'
                        event.Source.Tag = 'P';
                        PeriodicTable_PButtonPushed(app, event);
                end
                
                return

            % Modifier used (control + key press)
            elseif isequal(modifier, 'control')
                switch lower(key)
                    case 'n'
                        ClearButtonPushed(app, event);
                    case 'p'
                        PreferencesMenuSelected(app, event);
                    case 'q'
                        app.Console_text.Value = 'Closing the Combustion Toolbox';
                        UIFigureCloseRequest(app, event);
                    case 's'
                        xlsMenuSelected(app, event);
                    case 'x'
                        app.Console_text.Value = 'Taking snapshot...';
                        SnapshotMenuSelected(app, event);
                        app.Console_text.Value = 'Taking snapshot... OK!';
                end

                return
            
            % No modifier used (one key press)
            else
                switch lower(key)
                    case 'uparrow'
                        app.temp_index = app.temp_index - 1;
                    case 'downarrow'
                        app.temp_index = app.temp_index + 1;
                    case 'rightarrow'
                        [~, app.temp_index] = public_get_current_history(app);
                    case 'leftarrow'
                        app.temp_index = 1;
                    case 'return'
                        FLAG_ENTER = true;
                        gui_ConsoleValueChanged(app, event);
                    case 'f1'
                        DocumentationMenuSelected(app, event);
                    case 'f5'
                        CalculateButtonPushed(app, event);
                        return
                    otherwise
                        return
                end
            
            end
            % Print currhistory for that index
            if ~FLAG_ENTER
                try
                    app.Console.Value = char(app.current_history(app.temp_index));
                catch
                    app.Console.Value = 'Value out of current command history.';
                    [~, app.temp_index] = public_get_current_history(app);
                end
            end
        end

        % Menu selected function: NewMenu
        function NewMenuSelected(app, event)
            ClearButtonPushed(app, event);
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            delete(app)
        end

        % Value changing function: PR4
        function PR4ValueChanging(app, event)
            % Only for oblique/polar detonations
            if ~any(contains(app.ProblemType.Value, 'DET_OBLIQUE', 'IgnoreCase', true) | contains(app.ProblemType.Value, 'DET_POLAR', 'IgnoreCase', true))
                return
            end

            app.PP4.Value = '';
        end

        % Value changing function: PP4
        function PP4ValueChanging(app, event)
            % Only for oblique/polar detonations
            if ~any(contains(app.ProblemType.Value, 'DET_OBLIQUE', 'IgnoreCase', true) | contains(app.ProblemType.Value, 'DET_POLAR', 'IgnoreCase', true))
                return
            end

            app.PR4.Value = '';
        end

        % Value changing function: PR5
        function PR5ValueChanging(app, event)
            % Only for oblique/polar shocks
            if ~any(contains(app.ProblemType.Value, 'SHOCK_OBLIQUE', 'IgnoreCase', true) | contains(app.ProblemType.Value, 'SHOCK_POLAR', 'IgnoreCase', true))
                return
            end

            app.PP5.Value = '';
        end

        % Value changing function: PP5
        function PP5ValueChanging(app, event)
            % Only for oblique/polar shocks
            if ~any(contains(app.ProblemType.Value, 'SHOCK_OBLIQUE', 'IgnoreCase', true) | contains(app.ProblemType.Value, 'SHOCK_POLAR', 'IgnoreCase', true))
                return
            end

            app.PR5.Value = '';
        end

        % Value changed function: IonizedspeciesCheckBox
        function IonizedspeciesCheckBoxValueChanged(app, event)
            gui_update_ions(app);
        end

        % Callback function: Tree_variable_x
        function Tree_variable_xCheckedNodesChanged(app, event)
            % Do not allow selection of two or more variables in the x axis
            if isscalar(app.Tree_variable_x.CheckedNodes) || isempty(app.Tree_variable_x.CheckedNodes)
                return
            end

            if app.Tree_variable_x.CheckedNodes(1) == event.PreviousCheckedNodes
                app.Tree_variable_x.CheckedNodes(1) = [];
            else
                app.Tree_variable_x.CheckedNodes(2) = [];
            end
        end

        % Callback function: Tree_variable_y
        function Tree_variable_yCheckedNodesChanged(app, event)
            % If molar fraction or mass fraction is selected, do not allow
            % selection of two or more variables in the y axis
            
            if isscalar(app.Tree_variable_y.CheckedNodes) || isempty(app.Tree_variable_y.CheckedNodes)
                return
            end

            if sum(contains([app.Tree_variable_y.CheckedNodes.Text], {'Xi', 'Yi', 'molarFraction', 'massFraction'})) == 0
                return
            end
            
            if sum(contains([event.PreviousCheckedNodes.Text], {'Xi', 'Yi', 'molarFraction', 'massFraction'}))
                indexRemove = [];
                
                for i = 1:length(app.Tree_variable_y.CheckedNodes)
                    if any(strcmpi(app.Tree_variable_y.CheckedNodes(i).Text, {'Xi', 'Yi', 'molarFraction', 'massFraction'}))
                        indexRemove = [indexRemove, i];
                    end
                end

                app.Tree_variable_y.CheckedNodes(indexRemove) = [];
                return
            end

            if sum(contains([event.CheckedNodes.Text], {'Xi', 'Yi', 'molarFraction', 'massFraction'}))
                indexRemove = [];
                
                for i = 1:length(app.Tree_variable_y.CheckedNodes)
                    if ~strcmpi(app.Tree_variable_y.CheckedNodes(i).Text, {'Xi', 'Yi', 'molarFraction', 'massFraction'})
                        indexRemove = [indexRemove, i];
                    end
                end

                app.Tree_variable_y.CheckedNodes(indexRemove) = [];
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Color = [0.9098 0.9098 0.8902];
            app.UIFigure.Position = [650 100 639 768];
            app.UIFigure.Name = 'Combustion Toolbox';
            app.UIFigure.Icon = fullfile(pathToMLAPP, 'assets', 'logo_CT_noversion.png');
            app.UIFigure.Resize = 'off';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);
            app.UIFigure.KeyPressFcn = createCallbackFcn(app, @UIFigureKeyPress, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.UIFigure);
            app.FileMenu.Text = 'File';

            % Create NewMenu
            app.NewMenu = uimenu(app.FileMenu);
            app.NewMenu.MenuSelectedFcn = createCallbackFcn(app, @NewMenuSelected, true);
            app.NewMenu.Text = 'New';

            % Create OpenMenu
            app.OpenMenu = uimenu(app.FileMenu);
            app.OpenMenu.Enable = 'off';
            app.OpenMenu.Text = 'Open';

            % Create SaveasMenu
            app.SaveasMenu = uimenu(app.FileMenu);
            app.SaveasMenu.Separator = 'on';
            app.SaveasMenu.Text = 'Save as';

            % Create xlsMenu
            app.xlsMenu = uimenu(app.SaveasMenu);
            app.xlsMenu.MenuSelectedFcn = createCallbackFcn(app, @xlsMenuSelected, true);
            app.xlsMenu.Text = '*.xls';

            % Create matMenu
            app.matMenu = uimenu(app.SaveasMenu);
            app.matMenu.MenuSelectedFcn = createCallbackFcn(app, @matMenuSelected, true);
            app.matMenu.Text = '*.mat';

            % Create SnapshotMenu
            app.SnapshotMenu = uimenu(app.FileMenu);
            app.SnapshotMenu.MenuSelectedFcn = createCallbackFcn(app, @SnapshotMenuSelected, true);
            app.SnapshotMenu.Separator = 'on';
            app.SnapshotMenu.Text = 'Snapshot';

            % Create PreferencesMenu
            app.PreferencesMenu = uimenu(app.FileMenu);
            app.PreferencesMenu.MenuSelectedFcn = createCallbackFcn(app, @PreferencesMenuSelected, true);
            app.PreferencesMenu.Separator = 'on';
            app.PreferencesMenu.Text = 'Preferences';

            % Create HelpMenu
            app.HelpMenu = uimenu(app.UIFigure);
            app.HelpMenu.Text = 'Help';

            % Create DocumentationMenu
            app.DocumentationMenu = uimenu(app.HelpMenu);
            app.DocumentationMenu.MenuSelectedFcn = createCallbackFcn(app, @DocumentationMenuSelected, true);
            app.DocumentationMenu.Text = 'Documentation';

            % Create TutorialMenu
            app.TutorialMenu = uimenu(app.HelpMenu);
            app.TutorialMenu.Enable = 'off';
            app.TutorialMenu.Text = 'Tutorial';

            % Create ExamplesMenu
            app.ExamplesMenu = uimenu(app.HelpMenu);
            app.ExamplesMenu.Enable = 'off';
            app.ExamplesMenu.Text = 'Examples';

            % Create ValidationsMenu
            app.ValidationsMenu = uimenu(app.HelpMenu);
            app.ValidationsMenu.MenuSelectedFcn = createCallbackFcn(app, @ValidationsMenuSelected, true);
            app.ValidationsMenu.Text = 'Validations';

            % Create DiscussionsMenu
            app.DiscussionsMenu = uimenu(app.HelpMenu);
            app.DiscussionsMenu.MenuSelectedFcn = createCallbackFcn(app, @DiscussionsMenuSelected, true);
            app.DiscussionsMenu.Text = 'Discussions';

            % Create CheckforupdatesMenu
            app.CheckforupdatesMenu = uimenu(app.HelpMenu);
            app.CheckforupdatesMenu.MenuSelectedFcn = createCallbackFcn(app, @CheckforupdatesMenuSelected, true);
            app.CheckforupdatesMenu.Text = 'Check for updates';

            % Create LicenseMenu
            app.LicenseMenu = uimenu(app.HelpMenu);
            app.LicenseMenu.MenuSelectedFcn = createCallbackFcn(app, @LicenseMenuSelected, true);
            app.LicenseMenu.Separator = 'on';
            app.LicenseMenu.Text = 'License';

            % Create SendfeedbackMenu
            app.SendfeedbackMenu = uimenu(app.HelpMenu);
            app.SendfeedbackMenu.MenuSelectedFcn = createCallbackFcn(app, @SendfeedbackMenuSelected, true);
            app.SendfeedbackMenu.Separator = 'on';
            app.SendfeedbackMenu.Text = 'Send feedback';

            % Create AboutMenu
            app.AboutMenu = uimenu(app.HelpMenu);
            app.AboutMenu.MenuSelectedFcn = createCallbackFcn(app, @AboutMenuSelected, true);
            app.AboutMenu.Text = 'About...';

            % Create Tab_lateral_bar
            app.Tab_lateral_bar = uitabgroup(app.UIFigure);
            app.Tab_lateral_bar.AutoResizeChildren = 'off';
            app.Tab_lateral_bar.TabLocation = 'left';
            app.Tab_lateral_bar.Position = [0 -2 639 771];

            % Create SetupTab
            app.SetupTab = uitab(app.Tab_lateral_bar);
            app.SetupTab.AutoResizeChildren = 'off';
            app.SetupTab.Title = 'Setup';
            app.SetupTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.SetupTab);
            app.TabGroup2.AutoResizeChildren = 'off';
            app.TabGroup2.Position = [1 94 568 676];

            % Create InputsTab
            app.InputsTab = uitab(app.TabGroup2);
            app.InputsTab.AutoResizeChildren = 'off';
            app.InputsTab.Title = 'Inputs';
            app.InputsTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create DefinereactantsandspeciestobeconsideredPanel
            app.DefinereactantsandspeciestobeconsideredPanel = uipanel(app.InputsTab);
            app.DefinereactantsandspeciestobeconsideredPanel.AutoResizeChildren = 'off';
            app.DefinereactantsandspeciestobeconsideredPanel.BorderType = 'none';
            app.DefinereactantsandspeciestobeconsideredPanel.Title = 'Define reactants and species to be considered';
            app.DefinereactantsandspeciestobeconsideredPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.DefinereactantsandspeciestobeconsideredPanel.FontName = 'Arial';
            app.DefinereactantsandspeciestobeconsideredPanel.FontWeight = 'bold';
            app.DefinereactantsandspeciestobeconsideredPanel.Position = [9 329 558 314];

            % Create UITable_R
            app.UITable_R = uitable(app.DefinereactantsandspeciestobeconsideredPanel);
            app.UITable_R.ColumnName = {'Species'; 'Nº moles'; 'Mole fraction'; 'Type'; 'Temperature [K]'};
            app.UITable_R.RowName = {};
            app.UITable_R.ColumnSortable = true;
            app.UITable_R.ColumnEditable = [true true false true true];
            app.UITable_R.CellEditCallback = createCallbackFcn(app, @UITable_RCellEdit, true);
            app.UITable_R.Position = [3 1 547 171];

            % Create ReactantsPanel
            app.ReactantsPanel = uipanel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.ReactantsPanel.AutoResizeChildren = 'off';
            app.ReactantsPanel.BorderType = 'none';
            app.ReactantsPanel.Title = 'Reactants';
            app.ReactantsPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.ReactantsPanel.Position = [0 233 225 55];

            % Create Reactants
            app.Reactants = uidropdown(app.ReactantsPanel);
            app.Reactants.Items = {'', 'Air', 'Methane + Air', 'Ethane + Air', 'Propane + Air', 'Acetylene + Air', 'Ethylene + Air', 'Bencene + Air', 'Iso-octane + Air', 'Hydrogen + Air', 'Carbon monoxide + Air', 'Methanol + Air', 'Ethanol + Air', 'Natural Gas + Air', 'Syngas + Air', 'LH2 + LOX', 'RP1 + LOX'};
            app.Reactants.ItemsData = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30];
            app.Reactants.Editable = 'on';
            app.Reactants.ValueChangedFcn = createCallbackFcn(app, @ReactantsValueChanged, true);
            app.Reactants.BackgroundColor = [1 1 1];
            app.Reactants.Position = [4 6 222 25];
            app.Reactants.Value = 1;

            % Create ListofSpeciesPanel
            app.ListofSpeciesPanel = uipanel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.ListofSpeciesPanel.AutoResizeChildren = 'off';
            app.ListofSpeciesPanel.BorderType = 'none';
            app.ListofSpeciesPanel.Title = 'List of Species';
            app.ListofSpeciesPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.ListofSpeciesPanel.Position = [239 179 174 109];

            % Create listbox_Products
            app.listbox_Products = uilistbox(app.ListofSpeciesPanel);
            app.listbox_Products.Items = {};
            app.listbox_Products.Multiselect = 'on';
            app.listbox_Products.Position = [1 4 174 81];
            app.listbox_Products.Value = {};

            % Create OFEditFieldLabel
            app.OFEditFieldLabel = uilabel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.OFEditFieldLabel.HorizontalAlignment = 'right';
            app.OFEditFieldLabel.Position = [446 218 25 23];
            app.OFEditFieldLabel.Text = 'O/F';

            % Create edit_OF
            app.edit_OF = uieditfield(app.DefinereactantsandspeciestobeconsideredPanel, 'numeric');
            app.edit_OF.Position = [478 218 73 22];

            % Create edit_phi
            app.edit_phi = uieditfield(app.DefinereactantsandspeciestobeconsideredPanel, 'text');
            app.edit_phi.ValueChangedFcn = createCallbackFcn(app, @edit_phiValueChanged, true);
            app.edit_phi.HorizontalAlignment = 'center';
            app.edit_phi.Position = [478 185 73 22];
            app.edit_phi.Value = '-';

            % Create text_phi
            app.text_phi = uilabel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.text_phi.HorizontalAlignment = 'right';
            app.text_phi.Position = [445 183 26 25];
            app.text_phi.Text = 'Phi';

            % Create FuelEditFieldLabel
            app.FuelEditFieldLabel = uilabel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.FuelEditFieldLabel.HorizontalAlignment = 'right';
            app.FuelEditFieldLabel.Position = [411 252 60 22];
            app.FuelEditFieldLabel.Text = '% Fuel';

            % Create ProductsPanel
            app.ProductsPanel = uipanel(app.DefinereactantsandspeciestobeconsideredPanel);
            app.ProductsPanel.AutoResizeChildren = 'off';
            app.ProductsPanel.BorderType = 'none';
            app.ProductsPanel.Title = 'Products';
            app.ProductsPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.ProductsPanel.Position = [0 177 225 57];

            % Create Products
            app.Products = uidropdown(app.ProductsPanel);
            app.Products.Items = {'', 'Complete reaction', 'HC/O2/N2', 'HC/O2/N2 extended', 'HC/O2/N2 rich', 'HC/O2/N2 Propellants', 'Soot formation', 'Soot formation Extended', 'Hydrogen', 'Hydrogen (L)', 'Hydrogen ions', 'Dissociated air', 'Air ions'};
            app.Products.Editable = 'on';
            app.Products.ValueChangedFcn = createCallbackFcn(app, @ProductsValueChanged, true);
            app.Products.BackgroundColor = [1 1 1];
            app.Products.Position = [4 7 222 25];
            app.Products.Value = '';

            % Create PeriodicTable_P
            app.PeriodicTable_P = uibutton(app.ProductsPanel, 'push');
            app.PeriodicTable_P.ButtonPushedFcn = createCallbackFcn(app, @PeriodicTable_PButtonPushed, true);
            app.PeriodicTable_P.Tag = 'P';
            app.PeriodicTable_P.Icon = fullfile(pathToMLAPP, 'assets', 'logo_uielements.svg');
            app.PeriodicTable_P.IconAlignment = 'center';
            app.PeriodicTable_P.Tooltip = {''};
            app.PeriodicTable_P.Position = [190 40 36 18];
            app.PeriodicTable_P.Text = '';

            % Create edit_F
            app.edit_F = uieditfield(app.DefinereactantsandspeciestobeconsideredPanel, 'numeric');
            app.edit_F.Position = [478 252 73 22];

            % Create PeriodicTable_R
            app.PeriodicTable_R = uibutton(app.DefinereactantsandspeciestobeconsideredPanel, 'push');
            app.PeriodicTable_R.ButtonPushedFcn = createCallbackFcn(app, @PeriodicTable_RButtonPushed, true);
            app.PeriodicTable_R.Tag = 'R';
            app.PeriodicTable_R.Icon = fullfile(pathToMLAPP, 'assets', 'logo_uielements.svg');
            app.PeriodicTable_R.IconAlignment = 'center';
            app.PeriodicTable_R.Tooltip = {''};
            app.PeriodicTable_R.Position = [189 270 36 18];
            app.PeriodicTable_R.Text = '';

            % Create DefinestateofreactantsandproductsPanel
            app.DefinestateofreactantsandproductsPanel = uipanel(app.InputsTab);
            app.DefinestateofreactantsandproductsPanel.AutoResizeChildren = 'off';
            app.DefinestateofreactantsandproductsPanel.BorderType = 'none';
            app.DefinestateofreactantsandproductsPanel.Title = 'Define state of reactants and products';
            app.DefinestateofreactantsandproductsPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.DefinestateofreactantsandproductsPanel.FontName = 'Arial';
            app.DefinestateofreactantsandproductsPanel.FontWeight = 'bold';
            app.DefinestateofreactantsandproductsPanel.Position = [9 134 551 115];

            % Create PP2
            app.PP2 = uieditfield(app.DefinestateofreactantsandproductsPanel, 'text');
            app.PP2.ValueChangedFcn = createCallbackFcn(app, @PP2ValueChanged, true);
            app.PP2.HorizontalAlignment = 'center';
            app.PP2.FontWeight = 'bold';
            app.PP2.Position = [347 13 91 19];
            app.PP2.Value = '1';

            % Create PP1
            app.PP1 = uieditfield(app.DefinestateofreactantsandproductsPanel, 'text');
            app.PP1.HorizontalAlignment = 'center';
            app.PP1.FontWeight = 'bold';
            app.PP1.Position = [347 39 91 19];
            app.PP1.Value = '2500';

            % Create PR1
            app.PR1 = uieditfield(app.DefinestateofreactantsandproductsPanel, 'text');
            app.PR1.ValueChangedFcn = createCallbackFcn(app, @PR1ValueChanged, true);
            app.PR1.HorizontalAlignment = 'center';
            app.PR1.FontWeight = 'bold';
            app.PR1.Position = [60 39 91 19];
            app.PR1.Value = '300';

            % Create text_RP1
            app.text_RP1 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_RP1.HorizontalAlignment = 'center';
            app.text_RP1.Position = [161 39 176 19];
            app.text_RP1.Text = 'Temperature [K]';

            % Create text_R1
            app.text_R1 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_R1.HorizontalAlignment = 'center';
            app.text_R1.VerticalAlignment = 'top';
            app.text_R1.FontWeight = 'bold';
            app.text_R1.Position = [60 64 91 19];
            app.text_R1.Text = 'Reactants';

            % Create text_P1
            app.text_P1 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_P1.HorizontalAlignment = 'center';
            app.text_P1.VerticalAlignment = 'top';
            app.text_P1.FontWeight = 'bold';
            app.text_P1.Position = [310 64 165 19];
            app.text_P1.Text = 'Products';

            % Create PR2
            app.PR2 = uieditfield(app.DefinestateofreactantsandproductsPanel, 'text');
            app.PR2.ValueChangedFcn = createCallbackFcn(app, @PR2ValueChanged, true);
            app.PR2.HorizontalAlignment = 'center';
            app.PR2.FontWeight = 'bold';
            app.PR2.Position = [60 13 91 19];
            app.PR2.Value = '1';

            % Create text_RP2
            app.text_RP2 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_RP2.HorizontalAlignment = 'center';
            app.text_RP2.Position = [161 13 176 19];
            app.text_RP2.Text = 'Pressure [bar]';

            % Create text_RP1_2
            app.text_RP1_2 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_RP1_2.HorizontalAlignment = 'center';
            app.text_RP1_2.Visible = 'off';
            app.text_RP1_2.Position = [447 39 102 19];
            app.text_RP1_2.Text = 'Area ratio A_c/A_t';

            % Create text_RP2_2
            app.text_RP2_2 = uilabel(app.DefinestateofreactantsandproductsPanel);
            app.text_RP2_2.HorizontalAlignment = 'center';
            app.text_RP2_2.Visible = 'off';
            app.text_RP2_2.Position = [447 13 102 19];
            app.text_RP2_2.Text = 'Mass flux [kg/s]';

            % Create AdditionalconstraintsPanel
            app.AdditionalconstraintsPanel = uipanel(app.InputsTab);
            app.AdditionalconstraintsPanel.AutoResizeChildren = 'off';
            app.AdditionalconstraintsPanel.BorderType = 'none';
            app.AdditionalconstraintsPanel.Title = 'Additional constraints';
            app.AdditionalconstraintsPanel.Visible = 'off';
            app.AdditionalconstraintsPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.AdditionalconstraintsPanel.FontName = 'Arial';
            app.AdditionalconstraintsPanel.FontWeight = 'bold';
            app.AdditionalconstraintsPanel.Position = [9 1 551 133];

            % Create text_R2
            app.text_R2 = uilabel(app.AdditionalconstraintsPanel);
            app.text_R2.HorizontalAlignment = 'center';
            app.text_R2.VerticalAlignment = 'top';
            app.text_R2.FontWeight = 'bold';
            app.text_R2.Position = [60 84 91 19];
            app.text_R2.Text = 'Reactants';

            % Create text_RP
            app.text_RP = uilabel(app.AdditionalconstraintsPanel);
            app.text_RP.HorizontalAlignment = 'center';
            app.text_RP.VerticalAlignment = 'top';
            app.text_RP.FontWeight = 'bold';
            app.text_RP.Position = [154 84 190 19];
            app.text_RP.Text = 'Products';

            % Create PR3
            app.PR3 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PR3.ValueChangedFcn = createCallbackFcn(app, @PR3ValueChanged, true);
            app.PR3.ValueChangingFcn = createCallbackFcn(app, @PR3ValueChanging, true);
            app.PR3.HorizontalAlignment = 'center';
            app.PR3.FontWeight = 'bold';
            app.PR3.Position = [60 59 91 19];
            app.PR3.Value = '500';

            % Create text_RP3
            app.text_RP3 = uilabel(app.AdditionalconstraintsPanel);
            app.text_RP3.HorizontalAlignment = 'center';
            app.text_RP3.Position = [154 59 190 19];
            app.text_RP3.Text = 'Constant Enthalpy: hP = hR';

            % Create PR4
            app.PR4 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PR4.ValueChangedFcn = createCallbackFcn(app, @PR4ValueChanged, true);
            app.PR4.ValueChangingFcn = createCallbackFcn(app, @PR4ValueChanging, true);
            app.PR4.HorizontalAlignment = 'center';
            app.PR4.FontWeight = 'bold';
            app.PR4.Position = [60 33 91 19];
            app.PR4.Value = '1';

            % Create text_RP4
            app.text_RP4 = uilabel(app.AdditionalconstraintsPanel);
            app.text_RP4.HorizontalAlignment = 'center';
            app.text_RP4.Position = [154 33 190 19];
            app.text_RP4.Text = 'Pressure [bar]';

            % Create PP4
            app.PP4 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PP4.ValueChangingFcn = createCallbackFcn(app, @PP4ValueChanging, true);
            app.PP4.HorizontalAlignment = 'center';
            app.PP4.FontWeight = 'bold';
            app.PP4.Position = [347 33 91 19];
            app.PP4.Value = '1';

            % Create PP3
            app.PP3 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PP3.ValueChangingFcn = createCallbackFcn(app, @PP3ValueChanging, true);
            app.PP3.HorizontalAlignment = 'center';
            app.PP3.FontWeight = 'bold';
            app.PP3.Position = [347 59 91 19];
            app.PP3.Value = '2500';

            % Create text_P2
            app.text_P2 = uilabel(app.AdditionalconstraintsPanel);
            app.text_P2.HorizontalAlignment = 'center';
            app.text_P2.VerticalAlignment = 'top';
            app.text_P2.FontWeight = 'bold';
            app.text_P2.Position = [347 84 91 19];
            app.text_P2.Text = 'Products';

            % Create PR5
            app.PR5 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PR5.ValueChangingFcn = createCallbackFcn(app, @PR5ValueChanging, true);
            app.PR5.HorizontalAlignment = 'center';
            app.PR5.FontWeight = 'bold';
            app.PR5.Position = [60 7 91 19];
            app.PR5.Value = '1';

            % Create text_RP5
            app.text_RP5 = uilabel(app.AdditionalconstraintsPanel);
            app.text_RP5.HorizontalAlignment = 'center';
            app.text_RP5.Position = [172 7 155 19];
            app.text_RP5.Text = 'Wave/Deflection angle [deg]';

            % Create PP5
            app.PP5 = uieditfield(app.AdditionalconstraintsPanel, 'text');
            app.PP5.ValueChangingFcn = createCallbackFcn(app, @PP5ValueChanging, true);
            app.PP5.HorizontalAlignment = 'center';
            app.PP5.FontWeight = 'bold';
            app.PP5.Position = [347 7 91 19];
            app.PP5.Value = '1';

            % Create Clear
            app.Clear = uibutton(app.InputsTab, 'push');
            app.Clear.ButtonPushedFcn = createCallbackFcn(app, @ClearButtonPushed, true);
            app.Clear.Position = [484 14 70 25];
            app.Clear.Text = 'Clear';

            % Create Calculate
            app.Calculate = uibutton(app.InputsTab, 'push');
            app.Calculate.ButtonPushedFcn = createCallbackFcn(app, @CalculateButtonPushed, true);
            app.Calculate.Position = [484 46 70 25];
            app.Calculate.Text = 'Calculate';

            % Create SelectProblemTypePanel
            app.SelectProblemTypePanel = uipanel(app.InputsTab);
            app.SelectProblemTypePanel.AutoResizeChildren = 'off';
            app.SelectProblemTypePanel.BorderType = 'none';
            app.SelectProblemTypePanel.Title = 'Select Problem Type';
            app.SelectProblemTypePanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.SelectProblemTypePanel.FontName = 'Arial';
            app.SelectProblemTypePanel.FontWeight = 'bold';
            app.SelectProblemTypePanel.Position = [9 245 551 78];

            % Create ProblemType
            app.ProblemType = uidropdown(app.SelectProblemTypePanel);
            app.ProblemType.Items = {'TP:  Equilibrium composition at defined T and P', 'HP: Adiabatic T and composition at constant P', 'SP:  Isentropic compression/expansion to a specified P', 'TV: Equilibrium composition at defined T and constant V', 'EV: Adiabatic T and composition at constant V', 'SV: Isentropic compresion/expansion to a specified V', 'SHOCK_I: Planar incident shock wave', 'SHOCK_R: Planar reflected shock wave', 'SHOCK_OBLIQUE: Oblique incident shock wave', 'SHOCK_OBLIQUE_R: Oblique reflected shock wave', 'SHOCK_POLAR: Polar shocks', 'SHOCK_POLAR_R: Polar shocks (regular reflections)', 'DET: Chapman-Jouguet Detonation', 'DET_OVERDRIVEN: Overdriven detonation', 'DET_UNDERDRIVEN: Underdriven detonation', 'DET_R: Reflected Chapman-Jouguet detonation', 'DET_OVERDRIVEN_R: Overdriven reflected detonation', 'DET_UNDERDRIVEN_R: Underdriven reflected detonation', 'DET_OBLIQUE: Oblique incident detonation', 'DET_POLAR: Polar detonations', 'ROCKET: Rocket propellant  performance'};
            app.ProblemType.ItemsData = {'TP', 'HP', 'SP', 'TV', 'EV', 'SV', 'SHOCK_I', 'SHOCK_R', 'SHOCK_OBLIQUE', 'SHOCK_OBLIQUE_R', 'SHOCK_POLAR', 'SHOCK_POLAR_R', 'DET', 'DET_OVERDRIVEN', 'DET_UNDERDRIVEN', 'DET_R', 'DET_OVERDRIVEN_R', 'DET_UNDERDRIVEN_R', 'DET_OBLIQUE', 'DET_POLAR', 'ROCKET'};
            app.ProblemType.ValueChangedFcn = createCallbackFcn(app, @ProblemTypeValueChanged, true);
            app.ProblemType.Position = [10 18 403 25];
            app.ProblemType.Value = 'TP';

            % Create FrozenchemistryCheckBox
            app.FrozenchemistryCheckBox = uicheckbox(app.SelectProblemTypePanel);
            app.FrozenchemistryCheckBox.ValueChangedFcn = createCallbackFcn(app, @FrozenchemistryCheckBoxValueChanged, true);
            app.FrozenchemistryCheckBox.Text = 'Frozen chemistry';
            app.FrozenchemistryCheckBox.Position = [437 36 114 22];

            % Create IonizedspeciesCheckBox
            app.IonizedspeciesCheckBox = uicheckbox(app.SelectProblemTypePanel);
            app.IonizedspeciesCheckBox.ValueChangedFcn = createCallbackFcn(app, @IonizedspeciesCheckBoxValueChanged, true);
            app.IonizedspeciesCheckBox.Text = 'Ionized species';
            app.IonizedspeciesCheckBox.Position = [437 18 105 22];

            % Create FLAG_IAC
            app.FLAG_IAC = uicheckbox(app.SelectProblemTypePanel);
            app.FLAG_IAC.ValueChangedFcn = createCallbackFcn(app, @FLAG_IACValueChanged, true);
            app.FLAG_IAC.Visible = 'off';
            app.FLAG_IAC.Text = 'IAC model';
            app.FLAG_IAC.Position = [437 -1 78 22];
            app.FLAG_IAC.Value = true;

            % Create QuicksettingsTab
            app.QuicksettingsTab = uitab(app.TabGroup2);
            app.QuicksettingsTab.AutoResizeChildren = 'off';
            app.QuicksettingsTab.Title = 'Quick settings';
            app.QuicksettingsTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create TabGroup4
            app.TabGroup4 = uitabgroup(app.QuicksettingsTab);
            app.TabGroup4.AutoResizeChildren = 'off';
            app.TabGroup4.Position = [47 350 477 286];

            % Create SelectionofspeciesTab
            app.SelectionofspeciesTab = uitab(app.TabGroup4);
            app.SelectionofspeciesTab.AutoResizeChildren = 'off';
            app.SelectionofspeciesTab.Title = 'Selection of species';
            app.SelectionofspeciesTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create text_list_species
            app.text_list_species = uilabel(app.SelectionofspeciesTab);
            app.text_list_species.VerticalAlignment = 'top';
            app.text_list_species.Position = [10 223 176 22];
            app.text_list_species.Text = 'Search Species';

            % Create edit_seeker_1
            app.edit_seeker_1 = uieditfield(app.SelectionofspeciesTab, 'text');
            app.edit_seeker_1.ValueChangingFcn = createCallbackFcn(app, @edit_seeker_1ValueChanging, true);
            app.edit_seeker_1.Position = [10 190 176 26];

            % Create text_LS_DB
            app.text_LS_DB = uilabel(app.SelectionofspeciesTab);
            app.text_LS_DB.VerticalAlignment = 'top';
            app.text_LS_DB.Position = [10 160 176 22];
            app.text_LS_DB.Text = 'List of species DB';

            % Create listbox_LS_DB
            app.listbox_LS_DB = uilistbox(app.SelectionofspeciesTab);
            app.listbox_LS_DB.Items = {};
            app.listbox_LS_DB.Multiselect = 'on';
            app.listbox_LS_DB.Position = [10 22 176 131];
            app.listbox_LS_DB.Value = {};

            % Create AddButton1
            app.AddButton1 = uibutton(app.SelectionofspeciesTab, 'push');
            app.AddButton1.ButtonPushedFcn = createCallbackFcn(app, @AddButton1Pushed, true);
            app.AddButton1.Position = [196 130 70 23];
            app.AddButton1.Text = 'Add >>';

            % Create RemoveButton1
            app.RemoveButton1 = uibutton(app.SelectionofspeciesTab, 'push');
            app.RemoveButton1.ButtonPushedFcn = createCallbackFcn(app, @RemoveButton1Pushed, true);
            app.RemoveButton1.Position = [196 99 70 23];
            app.RemoveButton1.Text = 'Remove <<';

            % Create text_LS
            app.text_LS = uilabel(app.SelectionofspeciesTab);
            app.text_LS.VerticalAlignment = 'top';
            app.text_LS.Position = [276 160 129 22];
            app.text_LS.Text = 'List of species';

            % Create listbox_LS
            app.listbox_LS = uilistbox(app.SelectionofspeciesTab);
            app.listbox_LS.Items = {};
            app.listbox_LS.Multiselect = 'on';
            app.listbox_LS.Position = [276 22 176 131];
            app.listbox_LS.Value = {};

            % Create PeriodicTable_P_2
            app.PeriodicTable_P_2 = uibutton(app.SelectionofspeciesTab, 'push');
            app.PeriodicTable_P_2.ButtonPushedFcn = createCallbackFcn(app, @PeriodicTable_PButtonPushed, true);
            app.PeriodicTable_P_2.Tag = 'P';
            app.PeriodicTable_P_2.Icon = fullfile(pathToMLAPP, 'assets', 'logo_uielements.svg');
            app.PeriodicTable_P_2.IconAlignment = 'center';
            app.PeriodicTable_P_2.Tooltip = {''};
            app.PeriodicTable_P_2.Position = [150 227 36 18];
            app.PeriodicTable_P_2.Text = '';

            % Create DisplayspeciesTab
            app.DisplayspeciesTab = uitab(app.TabGroup4);
            app.DisplayspeciesTab.AutoResizeChildren = 'off';
            app.DisplayspeciesTab.Title = 'Display species';
            app.DisplayspeciesTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create text_list_species_2
            app.text_list_species_2 = uilabel(app.DisplayspeciesTab);
            app.text_list_species_2.VerticalAlignment = 'top';
            app.text_list_species_2.Position = [10 223 176 22];
            app.text_list_species_2.Text = 'Search Species';

            % Create edit_seeker_2
            app.edit_seeker_2 = uieditfield(app.DisplayspeciesTab, 'text');
            app.edit_seeker_2.ValueChangingFcn = createCallbackFcn(app, @edit_seeker_2ValueChanging, true);
            app.edit_seeker_2.Position = [10 190 176 26];

            % Create text_LS_2
            app.text_LS_2 = uilabel(app.DisplayspeciesTab);
            app.text_LS_2.VerticalAlignment = 'top';
            app.text_LS_2.Position = [10 160 176 22];
            app.text_LS_2.Text = 'List of species';

            % Create listbox_LS_2
            app.listbox_LS_2 = uilistbox(app.DisplayspeciesTab);
            app.listbox_LS_2.Items = {};
            app.listbox_LS_2.Multiselect = 'on';
            app.listbox_LS_2.Position = [10 22 176 131];
            app.listbox_LS_2.Value = {};

            % Create AddButton2
            app.AddButton2 = uibutton(app.DisplayspeciesTab, 'push');
            app.AddButton2.ButtonPushedFcn = createCallbackFcn(app, @AddButton2Pushed, true);
            app.AddButton2.Position = [196 130 70 23];
            app.AddButton2.Text = 'Add >>';

            % Create RemoveButton2
            app.RemoveButton2 = uibutton(app.DisplayspeciesTab, 'push');
            app.RemoveButton2.ButtonPushedFcn = createCallbackFcn(app, @RemoveButton2Pushed, true);
            app.RemoveButton2.Position = [196 99 70 23];
            app.RemoveButton2.Text = 'Remove <<';

            % Create listbox_LS_display
            app.listbox_LS_display = uilistbox(app.DisplayspeciesTab);
            app.listbox_LS_display.Items = {};
            app.listbox_LS_display.Multiselect = 'on';
            app.listbox_LS_display.Position = [276 22 176 131];
            app.listbox_LS_display.Value = {};

            % Create text_LS_display
            app.text_LS_display = uilabel(app.DisplayspeciesTab);
            app.text_LS_display.VerticalAlignment = 'top';
            app.text_LS_display.Position = [276 160 129 22];
            app.text_LS_display.Text = 'Display species';

            % Create TolerancesPanel
            app.TolerancesPanel = uipanel(app.QuicksettingsTab);
            app.TolerancesPanel.AutoResizeChildren = 'off';
            app.TolerancesPanel.BorderType = 'none';
            app.TolerancesPanel.Title = 'Tolerances';
            app.TolerancesPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.TolerancesPanel.FontWeight = 'bold';
            app.TolerancesPanel.Position = [49 174 229 166];

            % Create TraceoptionEditFieldLabel
            app.TraceoptionEditFieldLabel = uilabel(app.TolerancesPanel);
            app.TraceoptionEditFieldLabel.HorizontalAlignment = 'right';
            app.TraceoptionEditFieldLabel.Position = [43 110 91 22];
            app.TraceoptionEditFieldLabel.Text = 'Trace option';

            % Create TraceoptionEditField
            app.TraceoptionEditField = uieditfield(app.TolerancesPanel, 'numeric');
            app.TraceoptionEditField.Limits = [0 Inf];
            app.TraceoptionEditField.ValueDisplayFormat = '%.1e';
            app.TraceoptionEditField.ValueChangedFcn = createCallbackFcn(app, @TraceoptionEditFieldValueChanged, true);
            app.TraceoptionEditField.Position = [149 110 72 22];
            app.TraceoptionEditField.Value = 1e-14;

            % Create DisplaySpeciesEditFieldLabel
            app.DisplaySpeciesEditFieldLabel = uilabel(app.TolerancesPanel);
            app.DisplaySpeciesEditFieldLabel.HorizontalAlignment = 'right';
            app.DisplaySpeciesEditFieldLabel.Position = [43 13 91 22];
            app.DisplaySpeciesEditFieldLabel.Text = 'Display Species';

            % Create DisplaySpeciesEditField
            app.DisplaySpeciesEditField = uieditfield(app.TolerancesPanel, 'numeric');
            app.DisplaySpeciesEditField.ValueDisplayFormat = '%.1e';
            app.DisplaySpeciesEditField.ValueChangedFcn = createCallbackFcn(app, @DisplaySpeciesEditFieldValueChanged, true);
            app.DisplaySpeciesEditField.Position = [149 13 72 22];
            app.DisplaySpeciesEditField.Value = 1e-14;

            % Create RootFindingMethodEditFieldLabel
            app.RootFindingMethodEditFieldLabel = uilabel(app.TolerancesPanel);
            app.RootFindingMethodEditFieldLabel.HorizontalAlignment = 'right';
            app.RootFindingMethodEditFieldLabel.Position = [17 79 117 22];
            app.RootFindingMethodEditFieldLabel.Text = 'Root Finding Method';

            % Create RootFindingMethodEditField
            app.RootFindingMethodEditField = uieditfield(app.TolerancesPanel, 'numeric');
            app.RootFindingMethodEditField.Limits = [0 Inf];
            app.RootFindingMethodEditField.ValueDisplayFormat = '%.1e';
            app.RootFindingMethodEditField.ValueChangedFcn = createCallbackFcn(app, @RootFindingMethodEditFieldValueChanged, true);
            app.RootFindingMethodEditField.Position = [149 79 72 22];
            app.RootFindingMethodEditField.Value = 0.001;

            % Create ShocksandDetonationsEditFieldLabel
            app.ShocksandDetonationsEditFieldLabel = uilabel(app.TolerancesPanel);
            app.ShocksandDetonationsEditFieldLabel.HorizontalAlignment = 'right';
            app.ShocksandDetonationsEditFieldLabel.Position = [-2 47 136 22];
            app.ShocksandDetonationsEditFieldLabel.Text = 'Shocks and Detonations';

            % Create ShocksandDetonationsEditField
            app.ShocksandDetonationsEditField = uieditfield(app.TolerancesPanel, 'numeric');
            app.ShocksandDetonationsEditField.Limits = [0 Inf];
            app.ShocksandDetonationsEditField.ValueDisplayFormat = '%.1e';
            app.ShocksandDetonationsEditField.ValueChangedFcn = createCallbackFcn(app, @ShocksandDetonationsEditFieldValueChanged, true);
            app.ShocksandDetonationsEditField.Position = [149 47 72 22];
            app.ShocksandDetonationsEditField.Value = 5e-05;

            % Create TuningparametersPanel
            app.TuningparametersPanel = uipanel(app.QuicksettingsTab);
            app.TuningparametersPanel.AutoResizeChildren = 'off';
            app.TuningparametersPanel.BorderType = 'none';
            app.TuningparametersPanel.Title = 'Tuning parameters';
            app.TuningparametersPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.TuningparametersPanel.FontWeight = 'bold';
            app.TuningparametersPanel.Position = [308 143 216 197];

            % Create MaxiterationsRFMEditFieldLabel
            app.MaxiterationsRFMEditFieldLabel = uilabel(app.TuningparametersPanel);
            app.MaxiterationsRFMEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxiterationsRFMEditFieldLabel.Position = [17 141 110 22];
            app.MaxiterationsRFMEditFieldLabel.Text = 'Max iterations RFM';

            % Create MaxiterationsRFMEditField
            app.MaxiterationsRFMEditField = uieditfield(app.TuningparametersPanel, 'numeric');
            app.MaxiterationsRFMEditField.Limits = [1 Inf];
            app.MaxiterationsRFMEditField.ValueDisplayFormat = '%d';
            app.MaxiterationsRFMEditField.ValueChangedFcn = createCallbackFcn(app, @MaxiterationsRFMEditFieldValueChanged, true);
            app.MaxiterationsRFMEditField.Position = [142 141 72 22];
            app.MaxiterationsRFMEditField.Value = 30;

            % Create RFMT0_REditFieldLabel
            app.RFMT0_REditFieldLabel = uilabel(app.TuningparametersPanel);
            app.RFMT0_REditFieldLabel.HorizontalAlignment = 'right';
            app.RFMT0_REditFieldLabel.Position = [63 44 64 22];
            app.RFMT0_REditFieldLabel.Text = 'RFM T0_R';

            % Create RFMT0_REditField
            app.RFMT0_REditField = uieditfield(app.TuningparametersPanel, 'numeric');
            app.RFMT0_REditField.Limits = [0 Inf];
            app.RFMT0_REditField.ValueDisplayFormat = '%d K';
            app.RFMT0_REditField.ValueChangedFcn = createCallbackFcn(app, @RFMT0_REditFieldValueChanged, true);
            app.RFMT0_REditField.Position = [142 44 72 22];
            app.RFMT0_REditField.Value = 3000;

            % Create RFMT0_LEditFieldLabel
            app.RFMT0_LEditFieldLabel = uilabel(app.TuningparametersPanel);
            app.RFMT0_LEditFieldLabel.HorizontalAlignment = 'right';
            app.RFMT0_LEditFieldLabel.Position = [65 78 62 22];
            app.RFMT0_LEditFieldLabel.Text = 'RFM T0_L';

            % Create RFMT0_LEditField
            app.RFMT0_LEditField = uieditfield(app.TuningparametersPanel, 'numeric');
            app.RFMT0_LEditField.Limits = [0 Inf];
            app.RFMT0_LEditField.ValueDisplayFormat = '%d K';
            app.RFMT0_LEditField.ValueChangedFcn = createCallbackFcn(app, @RFMT0_LEditFieldValueChanged, true);
            app.RFMT0_LEditField.Position = [142 78 72 22];
            app.RFMT0_LEditField.Value = 1000;

            % Create RFMT0EditFieldLabel
            app.RFMT0EditFieldLabel = uilabel(app.TuningparametersPanel);
            app.RFMT0EditFieldLabel.HorizontalAlignment = 'right';
            app.RFMT0EditFieldLabel.Position = [78 9 49 22];
            app.RFMT0EditFieldLabel.Text = 'RFM T0';

            % Create RFMT0EditField
            app.RFMT0EditField = uieditfield(app.TuningparametersPanel, 'numeric');
            app.RFMT0EditField.Limits = [0 Inf];
            app.RFMT0EditField.ValueDisplayFormat = '%d K';
            app.RFMT0EditField.ValueChangedFcn = createCallbackFcn(app, @RFMT0EditFieldValueChanged, true);
            app.RFMT0EditField.Position = [142 9 72 22];
            app.RFMT0EditField.Value = 3000;

            % Create MaxiterationsSDEditFieldLabel
            app.MaxiterationsSDEditFieldLabel = uilabel(app.TuningparametersPanel);
            app.MaxiterationsSDEditFieldLabel.HorizontalAlignment = 'right';
            app.MaxiterationsSDEditFieldLabel.Position = [19 109 108 22];
            app.MaxiterationsSDEditFieldLabel.Text = 'Max iterations S&D';

            % Create MaxiterationsSDEditField
            app.MaxiterationsSDEditField = uieditfield(app.TuningparametersPanel, 'numeric');
            app.MaxiterationsSDEditField.Limits = [1 Inf];
            app.MaxiterationsSDEditField.ValueDisplayFormat = '%d';
            app.MaxiterationsSDEditField.ValueChangedFcn = createCallbackFcn(app, @MaxiterationsSDEditFieldValueChanged, true);
            app.MaxiterationsSDEditField.Position = [142 109 72 22];
            app.MaxiterationsSDEditField.Value = 50;

            % Create OthersPanel
            app.OthersPanel = uipanel(app.QuicksettingsTab);
            app.OthersPanel.AutoResizeChildren = 'off';
            app.OthersPanel.BorderType = 'none';
            app.OthersPanel.Title = 'Others';
            app.OthersPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.OthersPanel.FontWeight = 'bold';
            app.OthersPanel.Position = [48 13 100 131];

            % Create IdealAirCheckBox
            app.IdealAirCheckBox = uicheckbox(app.OthersPanel);
            app.IdealAirCheckBox.Text = 'Ideal Air';
            app.IdealAirCheckBox.Position = [1 76 100 22];

            % Create ReportPanel
            app.ReportPanel = uipanel(app.QuicksettingsTab);
            app.ReportPanel.AutoResizeChildren = 'off';
            app.ReportPanel.BorderType = 'none';
            app.ReportPanel.Title = 'Report';
            app.ReportPanel.BackgroundColor = [0.9098 0.9098 0.8902];
            app.ReportPanel.FontWeight = 'bold';
            app.ReportPanel.Position = [157 13 121 131];

            % Create TypeDropDownLabel
            app.TypeDropDownLabel = uilabel(app.ReportPanel);
            app.TypeDropDownLabel.HorizontalAlignment = 'right';
            app.TypeDropDownLabel.Position = [3 76 31 22];
            app.TypeDropDownLabel.Text = 'Type';

            % Create Report_type
            app.Report_type = uidropdown(app.ReportPanel);
            app.Report_type.Items = {'None', 'Auto'};
            app.Report_type.Position = [49 76 72 22];
            app.Report_type.Value = 'Auto';

            % Create PrintresultsCheckBox
            app.PrintresultsCheckBox = uicheckbox(app.ReportPanel);
            app.PrintresultsCheckBox.Text = 'Print results';
            app.PrintresultsCheckBox.Position = [6 40 115 22];

            % Create ResultsTab
            app.ResultsTab = uitab(app.Tab_lateral_bar);
            app.ResultsTab.AutoResizeChildren = 'off';
            app.ResultsTab.Title = 'Results';
            app.ResultsTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create Tab_results
            app.Tab_results = uitabgroup(app.ResultsTab);
            app.Tab_results.AutoResizeChildren = 'off';
            app.Tab_results.Position = [1 294 568 476];

            % Create ParametersTab
            app.ParametersTab = uitab(app.Tab_results);
            app.ParametersTab.AutoResizeChildren = 'off';
            app.ParametersTab.Title = 'Parameters';
            app.ParametersTab.BackgroundColor = [0.9098 0.9098 0.8902];
            app.ParametersTab.Scrollable = 'on';

            % Create Panel_parameters
            app.Panel_parameters = uipanel(app.ParametersTab);
            app.Panel_parameters.AutoResizeChildren = 'off';
            app.Panel_parameters.BorderType = 'none';
            app.Panel_parameters.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_parameters.Position = [92 6 475 441];

            % Create edit_phi3
            app.edit_phi3 = uieditfield(app.Panel_parameters, 'text');
            app.edit_phi3.Editable = 'off';
            app.edit_phi3.HorizontalAlignment = 'center';
            app.edit_phi3.Position = [162 407 64 19];
            app.edit_phi3.Value = '-';

            % Create text_phi_3
            app.text_phi_3 = uilabel(app.Panel_parameters);
            app.text_phi_3.HorizontalAlignment = 'center';
            app.text_phi_3.FontWeight = 'bold';
            app.text_phi_3.Position = [177 426 35 18];
            app.text_phi_3.Text = 'phi';

            % Create text_s
            app.text_s = uilabel(app.Panel_parameters);
            app.text_s.HorizontalAlignment = 'center';
            app.text_s.Position = [107 254 176 19];
            app.text_s.Text = 'Entropy [kJ/(kg K)]';

            % Create text_s_1
            app.text_s_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_s_1.ValueDisplayFormat = '%.4g';
            app.text_s_1.Editable = 'off';
            app.text_s_1.HorizontalAlignment = 'center';
            app.text_s_1.Position = [6 254 91 19];

            % Create text_s_2
            app.text_s_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_s_2.ValueDisplayFormat = '%.4g';
            app.text_s_2.Editable = 'off';
            app.text_s_2.HorizontalAlignment = 'center';
            app.text_s_2.Position = [293 254 91 19];

            % Create text_gamma
            app.text_gamma = uilabel(app.Panel_parameters);
            app.text_gamma.HorizontalAlignment = 'center';
            app.text_gamma.Position = [107 203 176 19];
            app.text_gamma.Text = 'Adiabatic index [-]';

            % Create text_gamma_1
            app.text_gamma_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_gamma_1.ValueDisplayFormat = '%.4g';
            app.text_gamma_1.Editable = 'off';
            app.text_gamma_1.HorizontalAlignment = 'center';
            app.text_gamma_1.Position = [6 203 91 19];

            % Create text_gamma_2
            app.text_gamma_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_gamma_2.ValueDisplayFormat = '%.4g';
            app.text_gamma_2.Editable = 'off';
            app.text_gamma_2.HorizontalAlignment = 'center';
            app.text_gamma_2.Position = [293 203 91 19];

            % Create text_M
            app.text_M = uilabel(app.Panel_parameters);
            app.text_M.HorizontalAlignment = 'center';
            app.text_M.Visible = 'off';
            app.text_M.Position = [107 100 176 19];
            app.text_M.Text = 'Mach number [-]';

            % Create text_M_2
            app.text_M_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_M_2.ValueDisplayFormat = '%.4g';
            app.text_M_2.Editable = 'off';
            app.text_M_2.HorizontalAlignment = 'center';
            app.text_M_2.Visible = 'off';
            app.text_M_2.Position = [293 100 91 19];

            % Create text_u
            app.text_u = uilabel(app.Panel_parameters);
            app.text_u.HorizontalAlignment = 'center';
            app.text_u.Visible = 'off';
            app.text_u.Position = [107 125 176 19];
            app.text_u.Text = 'Flow velocity [m/s]';

            % Create text_u_2
            app.text_u_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_u_2.ValueDisplayFormat = '%.4g';
            app.text_u_2.Editable = 'off';
            app.text_u_2.HorizontalAlignment = 'center';
            app.text_u_2.Visible = 'off';
            app.text_u_2.Position = [293 125 91 19];

            % Create text_e
            app.text_e = uilabel(app.Panel_parameters);
            app.text_e.HorizontalAlignment = 'center';
            app.text_e.Position = [107 280 176 19];
            app.text_e.Text = 'Internal energy [kJ/kg]';

            % Create text_e_1
            app.text_e_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_e_1.ValueDisplayFormat = '%.4g';
            app.text_e_1.Editable = 'off';
            app.text_e_1.HorizontalAlignment = 'center';
            app.text_e_1.Position = [6 280 91 19];

            % Create text_e_2
            app.text_e_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_e_2.ValueDisplayFormat = '%.4g';
            app.text_e_2.Editable = 'off';
            app.text_e_2.HorizontalAlignment = 'center';
            app.text_e_2.Position = [293 280 91 19];

            % Create text_sound
            app.text_sound = uilabel(app.Panel_parameters);
            app.text_sound.HorizontalAlignment = 'center';
            app.text_sound.Position = [107 151 176 19];
            app.text_sound.Text = 'Sound speed [m/s]';

            % Create text_sound_1
            app.text_sound_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_sound_1.ValueDisplayFormat = '%.4g';
            app.text_sound_1.Editable = 'off';
            app.text_sound_1.HorizontalAlignment = 'center';
            app.text_sound_1.Position = [6 151 91 19];

            % Create text_sound_2
            app.text_sound_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_sound_2.ValueDisplayFormat = '%.4g';
            app.text_sound_2.Editable = 'off';
            app.text_sound_2.HorizontalAlignment = 'center';
            app.text_sound_2.Position = [293 151 91 19];

            % Create text_W
            app.text_W = uilabel(app.Panel_parameters);
            app.text_W.HorizontalAlignment = 'center';
            app.text_W.Position = [107 174 176 22];
            app.text_W.Text = 'Molecular Weight [g/mol]';

            % Create text_W_1
            app.text_W_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_W_1.ValueDisplayFormat = '%.4g';
            app.text_W_1.Editable = 'off';
            app.text_W_1.HorizontalAlignment = 'center';
            app.text_W_1.Position = [6 177 91 19];

            % Create text_W_2
            app.text_W_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_W_2.ValueDisplayFormat = '%.4g';
            app.text_W_2.Editable = 'off';
            app.text_W_2.HorizontalAlignment = 'center';
            app.text_W_2.Position = [293 177 91 19];

            % Create text_cp
            app.text_cp = uilabel(app.Panel_parameters);
            app.text_cp.HorizontalAlignment = 'center';
            app.text_cp.Position = [107 228 176 19];
            app.text_cp.Text = 'Specific heat cp [kJ/(kg K)]';

            % Create text_cp_1
            app.text_cp_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_cp_1.ValueDisplayFormat = '%.4g';
            app.text_cp_1.Editable = 'off';
            app.text_cp_1.HorizontalAlignment = 'center';
            app.text_cp_1.Position = [6 228 91 19];

            % Create text_cp_2
            app.text_cp_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_cp_2.ValueDisplayFormat = '%.4g';
            app.text_cp_2.Editable = 'off';
            app.text_cp_2.HorizontalAlignment = 'center';
            app.text_cp_2.Position = [293 228 91 19];

            % Create text_h
            app.text_h = uilabel(app.Panel_parameters);
            app.text_h.HorizontalAlignment = 'center';
            app.text_h.Position = [107 305 176 19];
            app.text_h.Text = 'Enthalpy [kJ/kg]';

            % Create text_h_1
            app.text_h_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_h_1.ValueDisplayFormat = '%.4g';
            app.text_h_1.Editable = 'off';
            app.text_h_1.HorizontalAlignment = 'center';
            app.text_h_1.Position = [6 305 91 19];

            % Create text_h_2
            app.text_h_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_h_2.ValueDisplayFormat = '%.4g';
            app.text_h_2.Editable = 'off';
            app.text_h_2.HorizontalAlignment = 'center';
            app.text_h_2.Position = [293 305 91 19];

            % Create text_r
            app.text_r = uilabel(app.Panel_parameters);
            app.text_r.HorizontalAlignment = 'center';
            app.text_r.Position = [107 331 176 19];
            app.text_r.Text = 'Density [kg/m3]';

            % Create text_r_1
            app.text_r_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_r_1.ValueDisplayFormat = '%.4g';
            app.text_r_1.Editable = 'off';
            app.text_r_1.HorizontalAlignment = 'center';
            app.text_r_1.Position = [6 331 91 19];

            % Create text_r_2
            app.text_r_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_r_2.ValueDisplayFormat = '%.4g';
            app.text_r_2.Editable = 'off';
            app.text_r_2.HorizontalAlignment = 'center';
            app.text_r_2.Position = [293 331 91 19];

            % Create text_p
            app.text_p = uilabel(app.Panel_parameters);
            app.text_p.HorizontalAlignment = 'center';
            app.text_p.Position = [107 357 176 19];
            app.text_p.Text = 'Pressure [bar]';

            % Create text_p_2
            app.text_p_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_p_2.Limits = [0 Inf];
            app.text_p_2.ValueDisplayFormat = '%.4g';
            app.text_p_2.Editable = 'off';
            app.text_p_2.HorizontalAlignment = 'center';
            app.text_p_2.Position = [293 357 91 19];

            % Create text_Products
            app.text_Products = uilabel(app.Panel_parameters);
            app.text_Products.HorizontalAlignment = 'center';
            app.text_Products.FontWeight = 'bold';
            app.text_Products.Position = [293 408 91 19];
            app.text_Products.Text = 'Products';

            % Create text_Reactans
            app.text_Reactans = uilabel(app.Panel_parameters);
            app.text_Reactans.HorizontalAlignment = 'center';
            app.text_Reactans.FontWeight = 'bold';
            app.text_Reactans.Position = [6 405 91 22];
            app.text_Reactans.Text = 'Reactants';

            % Create text_T
            app.text_T = uilabel(app.Panel_parameters);
            app.text_T.HorizontalAlignment = 'center';
            app.text_T.Position = [107 383 176 19];
            app.text_T.Text = 'Temperature [K]';

            % Create text_T_1
            app.text_T_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_T_1.ValueDisplayFormat = '%.4g';
            app.text_T_1.Editable = 'off';
            app.text_T_1.HorizontalAlignment = 'center';
            app.text_T_1.Position = [6 383 91 19];

            % Create text_p_1
            app.text_p_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_p_1.ValueDisplayFormat = '%.4g';
            app.text_p_1.Editable = 'off';
            app.text_p_1.HorizontalAlignment = 'center';
            app.text_p_1.Position = [6 357 91 19];

            % Create text_T_2
            app.text_T_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_T_2.ValueDisplayFormat = '%.4g';
            app.text_T_2.Editable = 'off';
            app.text_T_2.HorizontalAlignment = 'center';
            app.text_T_2.Position = [293 383 91 19];

            % Create text_u_1
            app.text_u_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_u_1.ValueDisplayFormat = '%.4g';
            app.text_u_1.Editable = 'off';
            app.text_u_1.HorizontalAlignment = 'center';
            app.text_u_1.Visible = 'off';
            app.text_u_1.Position = [6 125 91 19];

            % Create text_M_1
            app.text_M_1 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_M_1.ValueDisplayFormat = '%.4g';
            app.text_M_1.Editable = 'off';
            app.text_M_1.HorizontalAlignment = 'center';
            app.text_M_1.Visible = 'off';
            app.text_M_1.Position = [6 100 91 19];

            % Create text_Aratio
            app.text_Aratio = uilabel(app.Panel_parameters);
            app.text_Aratio.HorizontalAlignment = 'center';
            app.text_Aratio.Visible = 'off';
            app.text_Aratio.Position = [129 75 136 19];
            app.text_Aratio.Text = 'Area ratio A/A_t [-]';

            % Create text_Aratio_2
            app.text_Aratio_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_Aratio_2.ValueDisplayFormat = '%.4g';
            app.text_Aratio_2.Editable = 'off';
            app.text_Aratio_2.HorizontalAlignment = 'center';
            app.text_Aratio_2.Visible = 'off';
            app.text_Aratio_2.Position = [293 75 91 19];

            % Create text_Cstar
            app.text_Cstar = uilabel(app.Panel_parameters);
            app.text_Cstar.HorizontalAlignment = 'center';
            app.text_Cstar.Visible = 'off';
            app.text_Cstar.Position = [137 51 120 19];
            app.text_Cstar.Text = 'Charac. velocity [m/s]';

            % Create text_Ivac
            app.text_Ivac = uilabel(app.Panel_parameters);
            app.text_Ivac.HorizontalAlignment = 'center';
            app.text_Ivac.Visible = 'off';
            app.text_Ivac.Position = [132 26 130 19];
            app.text_Ivac.Text = 'Specific impulse vac [s]';

            % Create text_Isp
            app.text_Isp = uilabel(app.Panel_parameters);
            app.text_Isp.HorizontalAlignment = 'center';
            app.text_Isp.Visible = 'off';
            app.text_Isp.Position = [143 1 108 19];
            app.text_Isp.Text = 'Specific impulse [s]';

            % Create Panel_extra_3
            app.Panel_extra_3 = uipanel(app.Panel_parameters);
            app.Panel_extra_3.AutoResizeChildren = 'off';
            app.Panel_extra_3.BorderType = 'none';
            app.Panel_extra_3.Visible = 'off';
            app.Panel_extra_3.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_extra_3.Position = [388 1 100 432];

            % Create text_s_3
            app.text_s_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_s_3.ValueDisplayFormat = '%.4g';
            app.text_s_3.Editable = 'off';
            app.text_s_3.HorizontalAlignment = 'center';
            app.text_s_3.Position = [5 254 91 19];

            % Create text_gamma_3
            app.text_gamma_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_gamma_3.ValueDisplayFormat = '%.4g';
            app.text_gamma_3.Editable = 'off';
            app.text_gamma_3.HorizontalAlignment = 'center';
            app.text_gamma_3.Position = [5 203 91 19];

            % Create text_M_3
            app.text_M_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_M_3.ValueDisplayFormat = '%.4g';
            app.text_M_3.Editable = 'off';
            app.text_M_3.HorizontalAlignment = 'center';
            app.text_M_3.Position = [5 100 91 19];

            % Create text_u_3
            app.text_u_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_u_3.ValueDisplayFormat = '%.4g';
            app.text_u_3.Editable = 'off';
            app.text_u_3.HorizontalAlignment = 'center';
            app.text_u_3.Position = [5 125 91 19];

            % Create text_e_3
            app.text_e_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_e_3.ValueDisplayFormat = '%.4g';
            app.text_e_3.Editable = 'off';
            app.text_e_3.HorizontalAlignment = 'center';
            app.text_e_3.Position = [5 280 91 19];

            % Create text_sound_3
            app.text_sound_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_sound_3.ValueDisplayFormat = '%.4g';
            app.text_sound_3.Editable = 'off';
            app.text_sound_3.HorizontalAlignment = 'center';
            app.text_sound_3.Position = [5 151 91 19];

            % Create text_W_3
            app.text_W_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_W_3.ValueDisplayFormat = '%.4g';
            app.text_W_3.Editable = 'off';
            app.text_W_3.HorizontalAlignment = 'center';
            app.text_W_3.Position = [5 177 91 19];

            % Create text_cp_3
            app.text_cp_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_cp_3.ValueDisplayFormat = '%.4g';
            app.text_cp_3.Editable = 'off';
            app.text_cp_3.HorizontalAlignment = 'center';
            app.text_cp_3.Position = [5 228 91 19];

            % Create text_h_3
            app.text_h_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_h_3.ValueDisplayFormat = '%.4g';
            app.text_h_3.Editable = 'off';
            app.text_h_3.HorizontalAlignment = 'center';
            app.text_h_3.Position = [5 305 91 19];

            % Create text_r_3
            app.text_r_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_r_3.ValueDisplayFormat = '%.4g';
            app.text_r_3.Editable = 'off';
            app.text_r_3.HorizontalAlignment = 'center';
            app.text_r_3.Position = [5 331 91 19];

            % Create text_p_3
            app.text_p_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_p_3.Limits = [0 Inf];
            app.text_p_3.ValueDisplayFormat = '%.4g';
            app.text_p_3.Editable = 'off';
            app.text_p_3.HorizontalAlignment = 'center';
            app.text_p_3.Position = [5 357 91 19];

            % Create text_Products_3
            app.text_Products_3 = uilabel(app.Panel_extra_3);
            app.text_Products_3.HorizontalAlignment = 'center';
            app.text_Products_3.VerticalAlignment = 'top';
            app.text_Products_3.FontWeight = 'bold';
            app.text_Products_3.Position = [5 405 91 22];
            app.text_Products_3.Text = 'Throat';

            % Create text_T_3
            app.text_T_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_T_3.ValueDisplayFormat = '%.4g';
            app.text_T_3.Editable = 'off';
            app.text_T_3.HorizontalAlignment = 'center';
            app.text_T_3.Position = [5 383 91 19];

            % Create text_Aratio_3
            app.text_Aratio_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_Aratio_3.ValueDisplayFormat = '%.4g';
            app.text_Aratio_3.Editable = 'off';
            app.text_Aratio_3.HorizontalAlignment = 'center';
            app.text_Aratio_3.Position = [5 75 91 19];

            % Create text_Cstar_3
            app.text_Cstar_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_Cstar_3.ValueDisplayFormat = '%.4g';
            app.text_Cstar_3.Editable = 'off';
            app.text_Cstar_3.HorizontalAlignment = 'center';
            app.text_Cstar_3.Position = [5 50 91 19];

            % Create text_Ivac_3
            app.text_Ivac_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_Ivac_3.ValueDisplayFormat = '%.4g';
            app.text_Ivac_3.Editable = 'off';
            app.text_Ivac_3.HorizontalAlignment = 'center';
            app.text_Ivac_3.Position = [5 26 91 19];

            % Create text_Isp_3
            app.text_Isp_3 = uieditfield(app.Panel_extra_3, 'numeric');
            app.text_Isp_3.ValueDisplayFormat = '%.4g';
            app.text_Isp_3.Editable = 'off';
            app.text_Isp_3.HorizontalAlignment = 'center';
            app.text_Isp_3.Position = [5 1 91 19];

            % Create Panel_extra_4
            app.Panel_extra_4 = uipanel(app.Panel_parameters);
            app.Panel_extra_4.AutoResizeChildren = 'off';
            app.Panel_extra_4.BorderType = 'none';
            app.Panel_extra_4.Visible = 'off';
            app.Panel_extra_4.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_extra_4.Position = [486 0 100 433];

            % Create text_s_4
            app.text_s_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_s_4.ValueDisplayFormat = '%.4g';
            app.text_s_4.Editable = 'off';
            app.text_s_4.HorizontalAlignment = 'center';
            app.text_s_4.Position = [3 255 91 19];

            % Create text_gamma_4
            app.text_gamma_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_gamma_4.ValueDisplayFormat = '%.4g';
            app.text_gamma_4.Editable = 'off';
            app.text_gamma_4.HorizontalAlignment = 'center';
            app.text_gamma_4.Position = [3 204 91 19];

            % Create text_M_4
            app.text_M_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_M_4.ValueDisplayFormat = '%.4g';
            app.text_M_4.Editable = 'off';
            app.text_M_4.HorizontalAlignment = 'center';
            app.text_M_4.Position = [3 101 91 19];

            % Create text_u_4
            app.text_u_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_u_4.ValueDisplayFormat = '%.4g';
            app.text_u_4.Editable = 'off';
            app.text_u_4.HorizontalAlignment = 'center';
            app.text_u_4.Position = [3 126 91 19];

            % Create text_e_4
            app.text_e_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_e_4.ValueDisplayFormat = '%.4g';
            app.text_e_4.Editable = 'off';
            app.text_e_4.HorizontalAlignment = 'center';
            app.text_e_4.Position = [3 281 91 19];

            % Create text_sound_4
            app.text_sound_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_sound_4.ValueDisplayFormat = '%.4g';
            app.text_sound_4.Editable = 'off';
            app.text_sound_4.HorizontalAlignment = 'center';
            app.text_sound_4.Position = [3 152 91 19];

            % Create text_W_4
            app.text_W_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_W_4.ValueDisplayFormat = '%.4g';
            app.text_W_4.Editable = 'off';
            app.text_W_4.HorizontalAlignment = 'center';
            app.text_W_4.Position = [3 178 91 19];

            % Create text_cp_4
            app.text_cp_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_cp_4.ValueDisplayFormat = '%.4g';
            app.text_cp_4.Editable = 'off';
            app.text_cp_4.HorizontalAlignment = 'center';
            app.text_cp_4.Position = [3 229 91 19];

            % Create text_h_4
            app.text_h_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_h_4.ValueDisplayFormat = '%.4g';
            app.text_h_4.Editable = 'off';
            app.text_h_4.HorizontalAlignment = 'center';
            app.text_h_4.Position = [3 306 91 19];

            % Create text_r_4
            app.text_r_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_r_4.ValueDisplayFormat = '%.4g';
            app.text_r_4.Editable = 'off';
            app.text_r_4.HorizontalAlignment = 'center';
            app.text_r_4.Position = [3 332 91 19];

            % Create text_p_4
            app.text_p_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_p_4.Limits = [0 Inf];
            app.text_p_4.ValueDisplayFormat = '%.4g';
            app.text_p_4.Editable = 'off';
            app.text_p_4.HorizontalAlignment = 'center';
            app.text_p_4.Position = [3 358 91 19];

            % Create text_Products_4
            app.text_Products_4 = uilabel(app.Panel_extra_4);
            app.text_Products_4.HorizontalAlignment = 'center';
            app.text_Products_4.VerticalAlignment = 'top';
            app.text_Products_4.FontWeight = 'bold';
            app.text_Products_4.Position = [3 406 91 22];
            app.text_Products_4.Text = 'Exit';

            % Create text_T_4
            app.text_T_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_T_4.ValueDisplayFormat = '%.4g';
            app.text_T_4.Editable = 'off';
            app.text_T_4.HorizontalAlignment = 'center';
            app.text_T_4.Position = [3 384 91 19];

            % Create text_Aratio_4
            app.text_Aratio_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_Aratio_4.ValueDisplayFormat = '%.4g';
            app.text_Aratio_4.Editable = 'off';
            app.text_Aratio_4.HorizontalAlignment = 'center';
            app.text_Aratio_4.Position = [3 76 91 19];

            % Create text_Cstar_4
            app.text_Cstar_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_Cstar_4.ValueDisplayFormat = '%.4g';
            app.text_Cstar_4.Editable = 'off';
            app.text_Cstar_4.HorizontalAlignment = 'center';
            app.text_Cstar_4.Position = [3 51 91 19];

            % Create text_Ivac_4
            app.text_Ivac_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_Ivac_4.ValueDisplayFormat = '%.4g';
            app.text_Ivac_4.Editable = 'off';
            app.text_Ivac_4.HorizontalAlignment = 'center';
            app.text_Ivac_4.Position = [3 27 91 19];

            % Create text_Isp_4
            app.text_Isp_4 = uieditfield(app.Panel_extra_4, 'numeric');
            app.text_Isp_4.ValueDisplayFormat = '%.4g';
            app.text_Isp_4.Editable = 'off';
            app.text_Isp_4.HorizontalAlignment = 'center';
            app.text_Isp_4.Position = [3 2 91 19];

            % Create Panel_extra_5
            app.Panel_extra_5 = uipanel(app.Panel_parameters);
            app.Panel_extra_5.AutoResizeChildren = 'off';
            app.Panel_extra_5.BorderType = 'none';
            app.Panel_extra_5.Visible = 'off';
            app.Panel_extra_5.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_extra_5.Position = [586 0 100 433];

            % Create text_s_5
            app.text_s_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_s_5.ValueDisplayFormat = '%.4g';
            app.text_s_5.Editable = 'off';
            app.text_s_5.HorizontalAlignment = 'center';
            app.text_s_5.Position = [5 255 91 19];

            % Create text_gamma_5
            app.text_gamma_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_gamma_5.ValueDisplayFormat = '%.4g';
            app.text_gamma_5.Editable = 'off';
            app.text_gamma_5.HorizontalAlignment = 'center';
            app.text_gamma_5.Position = [5 204 91 19];

            % Create text_M_5
            app.text_M_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_M_5.ValueDisplayFormat = '%.4g';
            app.text_M_5.Editable = 'off';
            app.text_M_5.HorizontalAlignment = 'center';
            app.text_M_5.Position = [5 101 91 19];

            % Create text_u_5
            app.text_u_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_u_5.ValueDisplayFormat = '%.4g';
            app.text_u_5.Editable = 'off';
            app.text_u_5.HorizontalAlignment = 'center';
            app.text_u_5.Position = [5 126 91 19];

            % Create text_e_5
            app.text_e_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_e_5.ValueDisplayFormat = '%.4g';
            app.text_e_5.Editable = 'off';
            app.text_e_5.HorizontalAlignment = 'center';
            app.text_e_5.Position = [5 281 91 19];

            % Create text_sound_5
            app.text_sound_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_sound_5.ValueDisplayFormat = '%.4g';
            app.text_sound_5.Editable = 'off';
            app.text_sound_5.HorizontalAlignment = 'center';
            app.text_sound_5.Position = [5 152 91 19];

            % Create text_W_5
            app.text_W_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_W_5.ValueDisplayFormat = '%.4g';
            app.text_W_5.Editable = 'off';
            app.text_W_5.HorizontalAlignment = 'center';
            app.text_W_5.Position = [5 178 91 19];

            % Create text_cp_5
            app.text_cp_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_cp_5.ValueDisplayFormat = '%.4g';
            app.text_cp_5.Editable = 'off';
            app.text_cp_5.HorizontalAlignment = 'center';
            app.text_cp_5.Position = [5 229 91 19];

            % Create text_h_5
            app.text_h_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_h_5.ValueDisplayFormat = '%.4g';
            app.text_h_5.Editable = 'off';
            app.text_h_5.HorizontalAlignment = 'center';
            app.text_h_5.Position = [5 306 91 19];

            % Create text_r_5
            app.text_r_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_r_5.ValueDisplayFormat = '%.4g';
            app.text_r_5.Editable = 'off';
            app.text_r_5.HorizontalAlignment = 'center';
            app.text_r_5.Position = [5 332 91 19];

            % Create text_p_5
            app.text_p_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_p_5.Limits = [0 Inf];
            app.text_p_5.ValueDisplayFormat = '%.4g';
            app.text_p_5.Editable = 'off';
            app.text_p_5.HorizontalAlignment = 'center';
            app.text_p_5.Position = [5 358 91 19];

            % Create text_Products_5
            app.text_Products_5 = uilabel(app.Panel_extra_5);
            app.text_Products_5.HorizontalAlignment = 'center';
            app.text_Products_5.VerticalAlignment = 'top';
            app.text_Products_5.FontWeight = 'bold';
            app.text_Products_5.Position = [5 406 91 22];
            app.text_Products_5.Text = 'Exit';

            % Create text_T_5
            app.text_T_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_T_5.ValueDisplayFormat = '%.4g';
            app.text_T_5.Editable = 'off';
            app.text_T_5.HorizontalAlignment = 'center';
            app.text_T_5.Position = [5 384 91 19];

            % Create text_Aratio_5
            app.text_Aratio_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_Aratio_5.ValueDisplayFormat = '%.4g';
            app.text_Aratio_5.Editable = 'off';
            app.text_Aratio_5.HorizontalAlignment = 'center';
            app.text_Aratio_5.Position = [5 76 91 19];

            % Create text_Cstar_5
            app.text_Cstar_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_Cstar_5.ValueDisplayFormat = '%.4g';
            app.text_Cstar_5.Editable = 'off';
            app.text_Cstar_5.HorizontalAlignment = 'center';
            app.text_Cstar_5.Position = [5 51 91 19];

            % Create text_Ivac_5
            app.text_Ivac_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_Ivac_5.ValueDisplayFormat = '%.4g';
            app.text_Ivac_5.Editable = 'off';
            app.text_Ivac_5.HorizontalAlignment = 'center';
            app.text_Ivac_5.Position = [5 27 91 19];

            % Create text_Isp_5
            app.text_Isp_5 = uieditfield(app.Panel_extra_5, 'numeric');
            app.text_Isp_5.ValueDisplayFormat = '%.4g';
            app.text_Isp_5.Editable = 'off';
            app.text_Isp_5.HorizontalAlignment = 'center';
            app.text_Isp_5.Position = [5 2 91 19];

            % Create text_error_problem
            app.text_error_problem = uieditfield(app.Panel_parameters, 'numeric');
            app.text_error_problem.ValueDisplayFormat = '%11.4e';
            app.text_error_problem.Editable = 'off';
            app.text_error_problem.Position = [6 26 91 19];

            % Create EpsilonmolesLabel_2
            app.EpsilonmolesLabel_2 = uilabel(app.Panel_parameters);
            app.EpsilonmolesLabel_2.HorizontalAlignment = 'center';
            app.EpsilonmolesLabel_2.Position = [-5 48 114 22];
            app.EpsilonmolesLabel_2.Text = 'Epsilon (method)';

            % Create text_Cstar_2
            app.text_Cstar_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_Cstar_2.ValueDisplayFormat = '%.4g';
            app.text_Cstar_2.Editable = 'off';
            app.text_Cstar_2.HorizontalAlignment = 'center';
            app.text_Cstar_2.Visible = 'off';
            app.text_Cstar_2.Position = [293 50 91 19];

            % Create text_Ivac_2
            app.text_Ivac_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_Ivac_2.ValueDisplayFormat = '%.4g';
            app.text_Ivac_2.Editable = 'off';
            app.text_Ivac_2.HorizontalAlignment = 'center';
            app.text_Ivac_2.Visible = 'off';
            app.text_Ivac_2.Position = [293 26 91 19];

            % Create text_Isp_2
            app.text_Isp_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_Isp_2.ValueDisplayFormat = '%.4g';
            app.text_Isp_2.Editable = 'off';
            app.text_Isp_2.HorizontalAlignment = 'center';
            app.text_Isp_2.Visible = 'off';
            app.text_Isp_2.Position = [293 1 91 19];

            % Create text_beta_min
            app.text_beta_min = uilabel(app.Panel_parameters);
            app.text_beta_min.HorizontalAlignment = 'center';
            app.text_beta_min.Visible = 'off';
            app.text_beta_min.Position = [129 75 136 19];
            app.text_beta_min.Text = 'Min wave angle [deg]';

            % Create text_beta
            app.text_beta = uilabel(app.Panel_parameters);
            app.text_beta.HorizontalAlignment = 'center';
            app.text_beta.Visible = 'off';
            app.text_beta.Position = [137 51 120 19];
            app.text_beta.Text = 'Wave angle [deg]';

            % Create text_theta
            app.text_theta = uilabel(app.Panel_parameters);
            app.text_theta.HorizontalAlignment = 'center';
            app.text_theta.Visible = 'off';
            app.text_theta.Position = [132 26 130 19];
            app.text_theta.Text = 'Deflection angle [deg]';

            % Create text_beta_min_2
            app.text_beta_min_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_beta_min_2.ValueDisplayFormat = '%.4g';
            app.text_beta_min_2.Editable = 'off';
            app.text_beta_min_2.HorizontalAlignment = 'center';
            app.text_beta_min_2.Visible = 'off';
            app.text_beta_min_2.Position = [293 75 91 19];

            % Create text_beta_2
            app.text_beta_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_beta_2.ValueDisplayFormat = '%.4g';
            app.text_beta_2.Editable = 'off';
            app.text_beta_2.HorizontalAlignment = 'center';
            app.text_beta_2.Visible = 'off';
            app.text_beta_2.Position = [293 50 91 19];

            % Create text_theta_2
            app.text_theta_2 = uieditfield(app.Panel_parameters, 'numeric');
            app.text_theta_2.ValueDisplayFormat = '%.4g';
            app.text_theta_2.Editable = 'off';
            app.text_theta_2.HorizontalAlignment = 'center';
            app.text_theta_2.Visible = 'off';
            app.text_theta_2.Position = [293 26 91 19];

            % Create MixturecompositionTab
            app.MixturecompositionTab = uitab(app.Tab_results);
            app.MixturecompositionTab.AutoResizeChildren = 'off';
            app.MixturecompositionTab.Title = 'Mixture composition';
            app.MixturecompositionTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create UITable_R2
            app.UITable_R2 = uitable(app.MixturecompositionTab);
            app.UITable_R2.ColumnName = {'Species'; 'Nº moles'; 'Mole fraction'};
            app.UITable_R2.RowName = {};
            app.UITable_R2.Position = [10 10 270 380];

            % Create UITable_P
            app.UITable_P = uitable(app.MixturecompositionTab);
            app.UITable_P.ColumnName = {'Species'; 'Nº moles'; 'Mole fraction'};
            app.UITable_P.RowName = {};
            app.UITable_P.Position = [290 10 270 380];

            % Create text_P3
            app.text_P3 = uilabel(app.MixturecompositionTab);
            app.text_P3.HorizontalAlignment = 'center';
            app.text_P3.VerticalAlignment = 'top';
            app.text_P3.FontWeight = 'bold';
            app.text_P3.Position = [379 402 91 19];
            app.text_P3.Text = 'Products';

            % Create text_R3
            app.text_R3 = uilabel(app.MixturecompositionTab);
            app.text_R3.HorizontalAlignment = 'center';
            app.text_R3.VerticalAlignment = 'top';
            app.text_R3.FontWeight = 'bold';
            app.text_R3.Position = [86 402 91 19];
            app.text_R3.Text = 'Reactants';

            % Create edit_phi2
            app.edit_phi2 = uieditfield(app.MixturecompositionTab, 'text');
            app.edit_phi2.Editable = 'off';
            app.edit_phi2.HorizontalAlignment = 'center';
            app.edit_phi2.Position = [252 402 64 19];
            app.edit_phi2.Value = '-';

            % Create text_phi_2
            app.text_phi_2 = uilabel(app.MixturecompositionTab);
            app.text_phi_2.HorizontalAlignment = 'center';
            app.text_phi_2.FontWeight = 'bold';
            app.text_phi_2.Position = [265 419 40 25];
            app.text_phi_2.Text = 'phi';

            % Create CustomFiguresTab
            app.CustomFiguresTab = uitab(app.Tab_results);
            app.CustomFiguresTab.AutoResizeChildren = 'off';
            app.CustomFiguresTab.Title = 'Custom Figures';
            app.CustomFiguresTab.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create UIAxes
            app.UIAxes = uiaxes(app.CustomFiguresTab);
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.LabelFontSizeMultiplier = 1;
            app.UIAxes.LineWidth = 1.2;
            app.UIAxes.Box = 'on';
            app.UIAxes.FontSize = 12;
            app.UIAxes.TitleFontSizeMultiplier = 1;
            app.UIAxes.Position = [8 10 552 260];

            % Create Tree_mixtures
            app.Tree_mixtures = uitree(app.CustomFiguresTab, 'checkbox');
            app.Tree_mixtures.Position = [12 317 150 125];

            % Create Mixtures
            app.Mixtures = uitreenode(app.Tree_mixtures);
            app.Mixtures.Text = 'Mixtures';

            % Create Tree_variable_x
            app.Tree_variable_x = uitree(app.CustomFiguresTab, 'checkbox');
            app.Tree_variable_x.Position = [181 317 180 125];

            % Create Variable_x
            app.Variable_x = uitreenode(app.Tree_variable_x);
            app.Variable_x.Text = 'Variable x';

            % Assign Checked Nodes
            app.Tree_variable_x.CheckedNodesChangedFcn = createCallbackFcn(app, @Tree_variable_xCheckedNodesChanged, true);

            % Create Tree_variable_y
            app.Tree_variable_y = uitree(app.CustomFiguresTab, 'checkbox');
            app.Tree_variable_y.Position = [379 317 180 125];

            % Create Variable_y
            app.Variable_y = uitreenode(app.Tree_variable_y);
            app.Variable_y.Text = 'Variable y';

            % Assign Checked Nodes
            app.Tree_variable_y.CheckedNodesChangedFcn = createCallbackFcn(app, @Tree_variable_yCheckedNodesChanged, true);

            % Create figure_plot
            app.figure_plot = uibutton(app.CustomFiguresTab, 'push');
            app.figure_plot.ButtonPushedFcn = createCallbackFcn(app, @figure_plotButtonPushed, true);
            app.figure_plot.Position = [411 279 70 25];
            app.figure_plot.Text = 'Plot';

            % Create figure_clear
            app.figure_clear = uibutton(app.CustomFiguresTab, 'push');
            app.figure_clear.ButtonPushedFcn = createCallbackFcn(app, @figure_clearButtonPushed, true);
            app.figure_clear.Position = [488 279 70 25];
            app.figure_clear.Text = 'Clear';

            % Create figure_size
            app.figure_size = uibutton(app.CustomFiguresTab, 'push');
            app.figure_size.ButtonPushedFcn = createCallbackFcn(app, @figure_sizeButtonPushed, true);
            app.figure_size.Position = [333 279 70 25];
            app.figure_size.Text = 'Maximize';

            % Create figure_settings
            app.figure_settings = uibutton(app.CustomFiguresTab, 'push');
            app.figure_settings.ButtonPushedFcn = createCallbackFcn(app, @figure_settingsButtonPushed, true);
            app.figure_settings.Position = [13 279 70 25];
            app.figure_settings.Text = 'Settings';

            % Create DefaultsettingsCheckBox
            app.DefaultsettingsCheckBox = uicheckbox(app.CustomFiguresTab);
            app.DefaultsettingsCheckBox.Text = 'Default settings';
            app.DefaultsettingsCheckBox.Position = [93 281 105 22];
            app.DefaultsettingsCheckBox.Value = true;

            % Create Tree
            app.Tree = uitree(app.ResultsTab);
            app.Tree.SelectionChangedFcn = createCallbackFcn(app, @TreeSelectionChanged, true);
            app.Tree.Position = [9 96 554 194];

            % Create Node_Results
            app.Node_Results = uitreenode(app.Tree);
            app.Node_Results.Text = 'Results';

            % Create Console
            app.Console = uitextarea(app.UIFigure);
            app.Console.ValueChangedFcn = createCallbackFcn(app, @ConsoleValueChanged, true);
            app.Console.Position = [98 6 500 20];

            % Create Label_Console
            app.Label_Console = uilabel(app.UIFigure);
            app.Label_Console.Position = [78 5 15 22];
            app.Label_Console.Text = '>>';

            % Create Lamp
            app.Lamp = uilamp(app.UIFigure);
            app.Lamp.Position = [606 4 23 23];
            app.Lamp.Color = [0.8 0.8 0.8];

            % Create Console_text
            app.Console_text = uitextarea(app.UIFigure);
            app.Console_text.Position = [78 30 554 57];
            app.Console_text.Value = {'Welcome to Combustion Toolbox vXX.XX.XX --- A MATLAB-GUI based open-source tool for solving gaseous combustion problems.'};

            % Create ContextMenu_CommandWindow
            app.ContextMenu_CommandWindow = uicontextmenu(app.UIFigure);

            % Create MaximizeMenu
            app.MaximizeMenu = uimenu(app.ContextMenu_CommandWindow);
            app.MaximizeMenu.MenuSelectedFcn = createCallbackFcn(app, @MaximizeMenuSelected, true);
            app.MaximizeMenu.Text = 'Maximize';

            % Create MinimizeMenu
            app.MinimizeMenu = uimenu(app.ContextMenu_CommandWindow);
            app.MinimizeMenu.MenuSelectedFcn = createCallbackFcn(app, @MinimizeMenuSelected, true);
            app.MinimizeMenu.Text = 'Minimize';

            % Create ClearCommandWindowMenu
            app.ClearCommandWindowMenu = uimenu(app.ContextMenu_CommandWindow);
            app.ClearCommandWindowMenu.MenuSelectedFcn = createCallbackFcn(app, @ClearCommandWindowMenuSelected, true);
            app.ClearCommandWindowMenu.Text = 'Clear Command Window';
            
            % Assign app.ContextMenu_CommandWindow
            app.Console.ContextMenu = app.ContextMenu_CommandWindow;
            app.Label_Console.ContextMenu = app.ContextMenu_CommandWindow;
            app.Lamp.ContextMenu = app.ContextMenu_CommandWindow;
            app.Console_text.ContextMenu = app.ContextMenu_CommandWindow;

            % Create ContextMenu_UITree
            app.ContextMenu_UITree = uicontextmenu(app.UIFigure);

            % Create RemoveMenu
            app.RemoveMenu = uimenu(app.ContextMenu_UITree);
            app.RemoveMenu.MenuSelectedFcn = createCallbackFcn(app, @RemoveMenuSelected, true);
            app.RemoveMenu.Text = 'Remove';

            % Create ContextMenu_UIAxes
            app.ContextMenu_UIAxes = uicontextmenu(app.UIFigure);

            % Create property_inspector_menu
            app.property_inspector_menu = uimenu(app.ContextMenu_UIAxes);
            app.property_inspector_menu.MenuSelectedFcn = createCallbackFcn(app, @property_inspector_menuSelected, true);
            app.property_inspector_menu.Text = 'Property inspector';

            % Create Menu
            app.Menu = uimenu(app.ContextMenu_UIAxes);
            
            % Assign app.ContextMenu_UIAxes
            app.UIAxes.ContextMenu = app.ContextMenu_UIAxes;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = combustion_toolbox

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIFigure)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.UIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
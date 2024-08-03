classdef uielements < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIElements             matlab.ui.Figure
        Panel_All              matlab.ui.container.Panel
        GridLayout             matlab.ui.container.GridLayout
        Panel_6                matlab.ui.container.Panel
        GridLayout4            matlab.ui.container.GridLayout
        condensed              matlab.ui.control.CheckBox
        burcat                 matlab.ui.control.CheckBox
        ionized                matlab.ui.control.CheckBox
        Panel_Close            matlab.ui.container.Panel
        Close                  matlab.ui.control.Button
        Panel_Export           matlab.ui.container.Panel
        Export                 matlab.ui.control.Button
        Database               matlab.ui.container.Panel
        GridLayout2            matlab.ui.container.GridLayout
        log                    matlab.ui.control.CheckBox
        type                   matlab.ui.control.DropDown
        temperature_range      matlab.ui.control.EditField
        property               matlab.ui.control.DropDown
        plot                   matlab.ui.control.Button
        RemoveButton1          matlab.ui.control.Button
        listbox_LS             matlab.ui.control.ListBox
        AddButton1             matlab.ui.control.Button
        listbox_LS_DB          matlab.ui.control.ListBox
        Copy                   matlab.ui.control.Button
        text_LS                matlab.ui.control.Label
        text_LS_DB             matlab.ui.control.Label
        edit_seeker            matlab.ui.control.EditField
        Panel_5                matlab.ui.container.Panel
        GridLayout3            matlab.ui.container.GridLayout
        text_comments          matlab.ui.control.TextArea
        CommentsTextAreaLabel  matlab.ui.control.Label
        text_ef                matlab.ui.control.NumericEditField
        text_ef_label          matlab.ui.control.Label
        text_hf                matlab.ui.control.NumericEditField
        text_hf_label          matlab.ui.control.Label
        text_W                 matlab.ui.control.NumericEditField
        text_W_label           matlab.ui.control.Label
        text_phase             matlab.ui.control.EditField
        text_phase_label       matlab.ui.control.Label
        text_codename          matlab.ui.control.EditField
        text_codename_label    matlab.ui.control.Label
        text_species           matlab.ui.control.EditField
        text_species_label     matlab.ui.control.Label
        text_list_species      matlab.ui.control.Label
        Legend                 matlab.ui.control.Image
        element_92             matlab.ui.control.Image
        element_91             matlab.ui.control.Image
        element_90             matlab.ui.control.Image
        element_89             matlab.ui.control.Image
        element_88             matlab.ui.control.Image
        element_87             matlab.ui.control.Image
        element_86             matlab.ui.control.Image
        element_85             matlab.ui.control.Image
        element_84             matlab.ui.control.Image
        element_83             matlab.ui.control.Image
        element_82             matlab.ui.control.Image
        element_81             matlab.ui.control.Image
        element_80             matlab.ui.control.Image
        element_79             matlab.ui.control.Image
        element_78             matlab.ui.control.Image
        element_77             matlab.ui.control.Image
        element_76             matlab.ui.control.Image
        element_75             matlab.ui.control.Image
        element_74             matlab.ui.control.Image
        element_73             matlab.ui.control.Image
        element_72             matlab.ui.control.Image
        element_57_71          matlab.ui.control.Image
        element_56             matlab.ui.control.Image
        element_55             matlab.ui.control.Image
        element_54             matlab.ui.control.Image
        element_53             matlab.ui.control.Image
        element_52             matlab.ui.control.Image
        element_51             matlab.ui.control.Image
        element_50             matlab.ui.control.Image
        element_49             matlab.ui.control.Image
        element_48             matlab.ui.control.Image
        element_47             matlab.ui.control.Image
        element_46             matlab.ui.control.Image
        element_45             matlab.ui.control.Image
        element_44             matlab.ui.control.Image
        element_43             matlab.ui.control.Image
        element_42             matlab.ui.control.Image
        element_41             matlab.ui.control.Image
        element_40             matlab.ui.control.Image
        element_39             matlab.ui.control.Image
        element_38             matlab.ui.control.Image
        element_37             matlab.ui.control.Image
        element_36             matlab.ui.control.Image
        element_35             matlab.ui.control.Image
        element_34             matlab.ui.control.Image
        element_33             matlab.ui.control.Image
        element_32             matlab.ui.control.Image
        element_31             matlab.ui.control.Image
        element_30             matlab.ui.control.Image
        element_29             matlab.ui.control.Image
        element_28             matlab.ui.control.Image
        element_27             matlab.ui.control.Image
        element_26             matlab.ui.control.Image
        element_25             matlab.ui.control.Image
        element_24             matlab.ui.control.Image
        element_23             matlab.ui.control.Image
        element_22             matlab.ui.control.Image
        element_21             matlab.ui.control.Image
        element_20             matlab.ui.control.Image
        element_19             matlab.ui.control.Image
        element_18             matlab.ui.control.Image
        element_17             matlab.ui.control.Image
        element_16             matlab.ui.control.Image
        element_15             matlab.ui.control.Image
        element_14             matlab.ui.control.Image
        element_13             matlab.ui.control.Image
        element_12             matlab.ui.control.Image
        element_11             matlab.ui.control.Image
        element_10             matlab.ui.control.Image
        element_9              matlab.ui.control.Image
        element_8              matlab.ui.control.Image
        element_7              matlab.ui.control.Image
        element_6              matlab.ui.control.Image
        element_5              matlab.ui.control.Image
        element_4              matlab.ui.control.Image
        element_3              matlab.ui.control.Image
        element_2              matlab.ui.control.Image
        element_1_3            matlab.ui.control.Image
        element_1_2            matlab.ui.control.Image
        element_1              matlab.ui.control.Image
        ContextMenu            matlab.ui.container.ContextMenu
        enableremoveMenu       matlab.ui.container.Menu
    end

    
    properties (Access = public)
        constants            % Constants class
        database             % Class inhereted from Database superclass
        chemicalSystem       % Chemical system class
        mixture              % Mixture class
        elements             % 
        plotConfig           % PlotConfig object
        indexElementsDatabase
        listElementsSelected % List of elements selected
        listElementsOmitted  % List of elements omitted
        FLAG_CALLER          % Flag app is called from other app
    end

    properties (Dependent)
        numElements
    end

    properties (Access = private)
        callerApp       % Handle to caller app
        callerAppTag    % Tag of caller app
        FLAG_LAST_CLICK % true == listbox_LS_DB, false = listbox_LS
    end
    
    methods (Access = private)
        
        function FLAG_CLICKED = changes_image(app, event)
            try
                image_source = split(event.Source.ImageSource, '\');
                image_element = image_source{end};

                if contains(image_element, 'clicked')
                    FLAG_CLICKED = false;
                    event.Source.ImageSource = strrep(image_element, '_clicked', '');
                else
                    FLAG_CLICKED = true;
                    event.Source.ImageSource = strrep(image_element, '.svg', '_clicked.svg');
                end

            catch
                FLAG_CLICKED = false;
            end

        end
        
        function plot_property(app)
            % Plot property against temperature
            
            % Import packages
            import combustiontoolbox.utils.display.*
            import combustiontoolbox.utils.extensions.brewermap

            % Get species
            if app.FLAG_LAST_CLICK
                species = app.listbox_LS_DB.Value;
            else
                species = app.listbox_LS.Value;
            end
            
            % Get number of species
            numSpecies = length(species);

            % Check that there are species selected
            if numSpecies == 0
                return
            end

            % Read temperature range
            T = gui_get_prop(app.temperature_range.Value);
            numTemperature = length(T);

            % Get function
            [funname_NASA, funname_CT, property_labelname] = set_inputs_thermo_validations(app.property.Value);
            
            if strcmpi(app.type.Value, 'CT')
                funname = funname_CT;
            else
                funname = funname_NASA;
            end

            % MATLAB OOP access to properties that contains a lot of data
            % is very very very ... slow.

            % tic
            % % Evaluate property for the temperature range
            % for i = numSpecies:-1:1
            %     for j = numTemperature:-1:1
            %         values(i, j) = funname(species{i}, temperature(j), app.database);
            %     end
            % 
            % end
            % toc
            
            % Initialization
            indexRemove = [];
            messages = [];
            
            % Copy DB local variables
            DB = app.database.species;

            % Evaluate property for the temperature range
            for i = numSpecies:-1:1
                try
                    for j = numTemperature:-1:1
                        values(i, j) = funname(species{i}, T(j), DB);
                    end

                catch
                    % Get index of the species that can not be evaluated
                    indexRemove = [indexRemove, i];
                    % Print message in the command window
                    message = ['%s can not be evaluated outside T = %.2f K.\n'...
                               'Species removed from the calculations.\n'];
                    message = sprintf(message, species{i}, DB.(species{i}).T);
                    % Update messages
                    messages = [messages, message, newline];
                    % Print
                    fprintf(message);
                end

                speciesLabel{i} = species2latex(species{i});
            end
            
            % Show warning if indexRemove
            if ~isempty(indexRemove)
                uialert(app.UIElements, messages, 'Warning', 'Icon', 'warning');
            end

            % Remove species that can not be evaluated
            species(indexRemove) = [];
            numSpecies = length(species);
            
            % Check that there are species that satisfied the requirements
            if numSpecies == 0
                return
            end
            
            % Get FLAG if the calculations imply extrapolation
            FLAG_EXTRAPOLATION_PRE = false(numSpecies, numTemperature);
            FLAG_EXTRAPOLATION_POST = false(numSpecies, numTemperature);
            for i = numSpecies:-1:1
                FLAG_EXTRAPOLATION_PRE(i, :) = T < DB.(species{i}).T(1);
                FLAG_EXTRAPOLATION_POST(i, :) = T > DB.(species{i}).T(end);
            end
            
            % Check figure scale
            FLAG_LOG = app.log.Value;
            if FLAG_LOG
                app.plotConfig.xscale = 'log';
                app.plotConfig.yscale = 'log';
            else
                app.plotConfig.xscale = 'linear';
                app.plotConfig.yscale = 'linear';
            end

            % Initialize figure
            ax = setFigure(app.plotConfig);
            
            % Define color palette
            color_palette = brewermap(numSpecies, app.plotConfig.colorpalette);

            % Plot
            dline = zeros(1, numSpecies);

            for i = 1:numSpecies
                try
                    [~, dline(i)] = plotFigure('T', T(~FLAG_EXTRAPOLATION_POST(i, :) & ~FLAG_EXTRAPOLATION_PRE(i, :)), property_labelname, values(i, ~FLAG_EXTRAPOLATION_POST(i, :) & ~FLAG_EXTRAPOLATION_PRE(i, :)), 'color', color_palette(i, :), 'linestyle', '-', 'ax', ax);
                catch
                    plotFigure('T', T(~FLAG_EXTRAPOLATION_POST(i, :) & ~FLAG_EXTRAPOLATION_PRE(i, :)), property_labelname, values(i, ~FLAG_EXTRAPOLATION_POST(i, :) & ~FLAG_EXTRAPOLATION_PRE(i, :)), 'color', color_palette(i, :), 'linestyle', '-', 'ax', ax);
                end
                
                try
                    [~, dline(i)] = plotFigure('T', T(FLAG_EXTRAPOLATION_PRE(i, :)), property_labelname, values(i, FLAG_EXTRAPOLATION_PRE(i, :)), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
                catch
                    plotFigure('T', T(FLAG_EXTRAPOLATION_PRE(i, :)), property_labelname, values(i, FLAG_EXTRAPOLATION_PRE(i, :)), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
                end

                try
                    [~, dline(i)] = plotFigure('T', T(FLAG_EXTRAPOLATION_POST(i, :)), property_labelname, values(i, FLAG_EXTRAPOLATION_POST(i, :)), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
                catch
                    plotFigure('T', T(FLAG_EXTRAPOLATION_POST(i, :)), property_labelname, values(i, FLAG_EXTRAPOLATION_POST(i, :)), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
                end
                
            end

            % Set legends
            if numSpecies < 2
                return
            end
            
            legend(ax, dline, speciesLabel, 'FontSize', app.plotConfig.fontsize - 4, 'Location', 'best', 'interpreter', 'latex')
        end

        
        function findSpecies(app)
            % Search species with the elements selected

            if app.numElements
                try
                    listSpecies = findProducts(app.chemicalSystem, app.listElementsSelected, 'ind', app.indexElementsDatabase);
                catch
                    listSpecies = {};
                end

                app.listbox_LS_DB.Items = listSpecies;
            end

        end

    end

    methods
        function value = get.numElements(app)
            value = length(app.listElementsSelected);
        end

    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, varargin)
            if nargin == 1
                app.FLAG_CALLER = false;
                % Get Constants
                app.constants = combustiontoolbox.common.Constants;
                % Get Nasa's database
                app.database = combustiontoolbox.databases.NasaDatabase();
                % Initialize chemical system
                app.chemicalSystem = combustiontoolbox.core.ChemicalSystem(app.database);
                % Initialize mixture
                app.mixture = combustiontoolbox.core.Mixture(app.chemicalSystem);
                % Initialize plotConfig
                app.plotConfig = combustiontoolbox.utils.display.PlotConfig();
            else
                app.FLAG_CALLER = true;
                app.callerApp = varargin{1};
                app.callerAppTag = varargin{2}.Source.Tag;
                app.constants = app.callerApp.constants;
                app.database = app.callerApp.database;
                app.chemicalSystem = app.callerApp.chemicalSystem;
                app.mixture = app.callerApp.mixture;
                app.plotConfig = app.callerApp.plotConfig;
            end
            
            % Initialize elements
            app.elements = combustiontoolbox.core.Elements();
            
            % Initialization
            app.listElementsOmitted = [];
            
            % Remove incompatible species
            if isfield(app.database.species, 'Air')
                app.database.species = rmfield(app.database.species, 'Air');
            end

            % Get element indeces of each species contained in the database
            app.indexElementsDatabase = getIndexElements(app.database.listSpecies, app.database.species, app.elements.listElements, 5);

            % Default size for plots
            % app.plotConfig.innerposition = [0.2, 0.15, 0.7, 0.7];
            % app.plotConfig.outerposition = [0.2, 0.15, 0.7, 0.7];

            % Set default values
            app.burcat.Value = app.chemicalSystem.FLAG_BURCAT;
            app.ionized.Value = app.chemicalSystem.FLAG_ION;
            app.condensed.Value = app.chemicalSystem.FLAG_CONDENSED;
        end

        % Image clicked function: element_1, element_10, element_11, 
        % ...and 77 other components
        function element_1ImageClicked(app, event)
            % Change image

            % Import pakcages
            import combustiontoolbox.utils.findIndex

            FLAG_CLICKED = changes_image(app, event);
            if FLAG_CLICKED
                % Add element to the list of elements
                if app.numElements
                    app.listElementsSelected{end + 1} = event.Source.Tag;
                else
                    app.listElementsSelected = {event.Source.Tag};
                end

            else
                % Remove element of the list of elements
                index_remove = findIndex(app.listElementsSelected, event.Source.Tag);
                app.listElementsSelected(index_remove) = [];
                if ~app.numElements
                    app.listbox_LS_DB.Items = {};
                end

            end

            % Search species with the elements selected
            findSpecies(app);
        end

        % Button pushed function: AddButton1
        function AddButton1Pushed(app, event)
            app.listbox_LS.Items = gui_value2list(app, app.listbox_LS_DB.Value, app.listbox_LS.Items, 'add');
        end

        % Button pushed function: RemoveButton1
        function RemoveButton1Pushed(app, event)
            app.listbox_LS.Items = gui_value2list(app, app.listbox_LS.Value, app.listbox_LS.Items, 'remove');
        end

        % Value changed function: listbox_LS, listbox_LS_DB
        function listbox_LSValueChanged(app, event)
            species = event.Source.Value{1};
            % Get values
            phase = app.database.species.(species).phase;
            W = app.database.species.(species).W; % [kg/mol]
            % Update GUI
            app.text_species.Value = app.database.species.(species).fullname;
            app.text_codename.Value = app.database.species.(species).name;
            app.text_comments.Value = strtrim(app.database.species.(species).comments);

            if phase
                app.text_phase.Value = 'condensed';
            else
                app.text_phase.Value = 'gas';
            end

            app.text_W.Value = W * 1e3; % [g/mol]
            app.text_hf.Value = app.database.species.(species).hf * W * 1e-3; % [kJ]
            app.text_ef.Value = app.database.species.(species).ef * W * 1e-3; % [kJ]
        end

        % Button pushed function: Export
        function ExportButtonPushed(app, event)
            if app.FLAG_CALLER
                % Passing data between apps
                switch upper(app.callerAppTag)
                    case 'R'
                        % To do it compatible with the current routines in the main app
                        app.callerApp.Reactants.Value = ' ';
                        % Update main GUI
                        for i = 1:length(app.listbox_LS.Items)
                            temp.Value = app.listbox_LS.Items{i};
                            public_ReactantsValueChanged(app.callerApp, temp)
                        end
                        
                    case 'P'
                        % Update main GUI
                        app.callerApp.listbox_Products.Items = unique([app.callerApp.listbox_Products.Items, app.listbox_LS.Items], 'stable');
                        public_ProductsValueChanged(app.callerApp)
                    otherwise
                        CopyButtonPushed(app, event);
                end
            else
                % Export as a txt and to the clipboard
                CopyButtonPushed(app, event);
            end
        end

        % Button pushed function: Copy
        function CopyButtonPushed(app, event)
            % Generate cell
            LS_copy = '{';
            for i = 1:length(app.listbox_LS.Items)
                if i > 1
                    LS_copy = [LS_copy, ', '];
                end

                LS_copy = [LS_copy, '''', app.listbox_LS.Items{i}, ''''];
            end

            LS_copy = [LS_copy, '}'];
            % Copy to clipboard
            clipboard('copy', LS_copy);
        end

        % Menu selected function: enableremoveMenu
        function enableremoveMenuSelected(app, event)
            % Get current object

            % Import pakcages
            import combustiontoolbox.utils.findIndex
            
            current_fig = gcbf;
            current_obj = current_fig.CurrentObject;
            % Check state to enable/remove the element
            if strcmpi(current_obj.Enable, 'on')
                % Omit element to be searchable
                current_obj.Enable = 'off';
                % Add element to listElementsOmitted
                app.listElementsOmitted{end+1} = current_obj.Tag; 
            else
                % Allow element to be searchable
                current_obj.Enable = 'on';
                % Remove element from listElementsOmitted
                ind = findIndex(app.listElementsOmitted, current_obj.Tag);
                app.listElementsOmitted(ind) = [];
            end
            
            % Update GUI
            element_1ImageClicked(app, event);
        end

        % Callback function: Close, UIElements
        function UIElementsCloseRequest(app, event)
            delete(app)
        end

        % Value changing function: edit_seeker
        function edit_seekerValueChanging(app, event)
            seek_value = gui_seeker_value(app, event, app.S.LS_DB);
            % Update Listbox (inputs)
            if ~isempty(seek_value)
                app.listbox_LS_DB.Items = seek_value;
            else
                app.listbox_LS_DB.Items = {};
            end
        end

        % Button pushed function: plot
        function plotButtonPushed(app, event)
            % Plot property against temperature
            plot_property(app);
        end

        % Clicked callback: listbox_LS_DB
        function listbox_LS_DBClicked(app, event)
            % Get last click | true == listbox_LS_DB, false = listbox_LS
            app.FLAG_LAST_CLICK = true;
        end

        % Clicked callback: listbox_LS
        function listbox_LSClicked(app, event)
            % Get last click | true == listbox_LS_DB, false = listbox_LS
            app.FLAG_LAST_CLICK = false;
        end

        % Value changed function: ionized
        function ionizedValueChanged(app, event)
            % Consider ionized species
            app.chemicalSystem.FLAG_ION = app.ionized.Value;
            % Search species
            findSpecies(app);
        end

        % Value changed function: condensed
        function condensedValueChanged(app, event)
            % Consider condensed species
            app.chemicalSystem.FLAG_CONDENSED = app.condensed.Value;
            % Search species
            findSpecies(app);
        end

        % Value changed function: burcat
        function burcatValueChanged(app, event)
            % Consider species from the Third Millenium database (Burcat)
            app.chemicalSystem.FLAG_BURCAT = app.burcat.Value;
            % Search species
            findSpecies(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIElements and hide until all components are created
            app.UIElements = uifigure('Visible', 'off');
            app.UIElements.Color = [0.9098 0.9098 0.8902];
            app.UIElements.Position = [250 250 1434 558];
            app.UIElements.Name = 'uielements';
            app.UIElements.Icon = 'logo_uielements.png';
            app.UIElements.CloseRequestFcn = createCallbackFcn(app, @UIElementsCloseRequest, true);
            app.UIElements.Scrollable = 'on';

            % Create Panel_All
            app.Panel_All = uipanel(app.UIElements);
            app.Panel_All.BorderType = 'none';
            app.Panel_All.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_All.Position = [1 1 1434 558];

            % Create GridLayout
            app.GridLayout = uigridlayout(app.Panel_All);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            app.GridLayout.Padding = [9.04736328125 8.19999186197917 9.04736328125 8.19999186197917];
            app.GridLayout.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create element_1
            app.element_1 = uiimage(app.GridLayout);
            app.element_1.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_1.Tag = 'H';
            app.element_1.Layout.Row = 1;
            app.element_1.Layout.Column = 1;
            app.element_1.ImageSource = '1.svg';

            % Create element_1_2
            app.element_1_2 = uiimage(app.GridLayout);
            app.element_1_2.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_1_2.Tag = 'D';
            app.element_1_2.Layout.Row = 1;
            app.element_1_2.Layout.Column = 2;
            app.element_1_2.ImageSource = '1_2.svg';

            % Create element_1_3
            app.element_1_3 = uiimage(app.GridLayout);
            app.element_1_3.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_1_3.Tag = 'T_M';
            app.element_1_3.Layout.Row = 1;
            app.element_1_3.Layout.Column = 3;
            app.element_1_3.ImageSource = '1_3.svg';

            % Create element_2
            app.element_2 = uiimage(app.GridLayout);
            app.element_2.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_2.Tag = 'He';
            app.element_2.Layout.Row = 1;
            app.element_2.Layout.Column = 18;
            app.element_2.ImageSource = '2.svg';

            % Create element_3
            app.element_3 = uiimage(app.GridLayout);
            app.element_3.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_3.Tag = 'Li';
            app.element_3.Layout.Row = 2;
            app.element_3.Layout.Column = 1;
            app.element_3.ImageSource = '3.svg';

            % Create element_4
            app.element_4 = uiimage(app.GridLayout);
            app.element_4.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_4.Tag = 'Be';
            app.element_4.Layout.Row = 2;
            app.element_4.Layout.Column = 2;
            app.element_4.ImageSource = '4.svg';

            % Create element_5
            app.element_5 = uiimage(app.GridLayout);
            app.element_5.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_5.Tag = 'B';
            app.element_5.Layout.Row = 2;
            app.element_5.Layout.Column = 13;
            app.element_5.ImageSource = '5.svg';

            % Create element_6
            app.element_6 = uiimage(app.GridLayout);
            app.element_6.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_6.Tag = 'C';
            app.element_6.Layout.Row = 2;
            app.element_6.Layout.Column = 14;
            app.element_6.ImageSource = '6.svg';

            % Create element_7
            app.element_7 = uiimage(app.GridLayout);
            app.element_7.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_7.Tag = 'N';
            app.element_7.Layout.Row = 2;
            app.element_7.Layout.Column = 15;
            app.element_7.ImageSource = '7.svg';

            % Create element_8
            app.element_8 = uiimage(app.GridLayout);
            app.element_8.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_8.Tag = 'O';
            app.element_8.Layout.Row = 2;
            app.element_8.Layout.Column = 16;
            app.element_8.ImageSource = '8.svg';

            % Create element_9
            app.element_9 = uiimage(app.GridLayout);
            app.element_9.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_9.Tag = 'F';
            app.element_9.Layout.Row = 2;
            app.element_9.Layout.Column = 17;
            app.element_9.ImageSource = '9.svg';

            % Create element_10
            app.element_10 = uiimage(app.GridLayout);
            app.element_10.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_10.Tag = 'Ne';
            app.element_10.Layout.Row = 2;
            app.element_10.Layout.Column = 18;
            app.element_10.ImageSource = '10.svg';

            % Create element_11
            app.element_11 = uiimage(app.GridLayout);
            app.element_11.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_11.Tag = 'Na';
            app.element_11.Layout.Row = 3;
            app.element_11.Layout.Column = 1;
            app.element_11.ImageSource = '11.svg';

            % Create element_12
            app.element_12 = uiimage(app.GridLayout);
            app.element_12.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_12.Tag = 'Mg';
            app.element_12.Layout.Row = 3;
            app.element_12.Layout.Column = 2;
            app.element_12.ImageSource = '12.svg';

            % Create element_13
            app.element_13 = uiimage(app.GridLayout);
            app.element_13.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_13.Tag = 'AL';
            app.element_13.Layout.Row = 3;
            app.element_13.Layout.Column = 13;
            app.element_13.ImageSource = '13.svg';

            % Create element_14
            app.element_14 = uiimage(app.GridLayout);
            app.element_14.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_14.Tag = 'Si';
            app.element_14.Layout.Row = 3;
            app.element_14.Layout.Column = 14;
            app.element_14.ImageSource = '14.svg';

            % Create element_15
            app.element_15 = uiimage(app.GridLayout);
            app.element_15.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_15.Tag = 'P';
            app.element_15.Layout.Row = 3;
            app.element_15.Layout.Column = 15;
            app.element_15.ImageSource = '15.svg';

            % Create element_16
            app.element_16 = uiimage(app.GridLayout);
            app.element_16.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_16.Tag = 'S';
            app.element_16.Layout.Row = 3;
            app.element_16.Layout.Column = 16;
            app.element_16.ImageSource = '16.svg';

            % Create element_17
            app.element_17 = uiimage(app.GridLayout);
            app.element_17.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_17.Tag = 'CL';
            app.element_17.Layout.Row = 3;
            app.element_17.Layout.Column = 17;
            app.element_17.ImageSource = '17.svg';

            % Create element_18
            app.element_18 = uiimage(app.GridLayout);
            app.element_18.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_18.Tag = 'Ar';
            app.element_18.Layout.Row = 3;
            app.element_18.Layout.Column = 18;
            app.element_18.ImageSource = '18.svg';

            % Create element_19
            app.element_19 = uiimage(app.GridLayout);
            app.element_19.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_19.Tag = 'K';
            app.element_19.Layout.Row = 4;
            app.element_19.Layout.Column = 1;
            app.element_19.ImageSource = '19.svg';

            % Create element_20
            app.element_20 = uiimage(app.GridLayout);
            app.element_20.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_20.Tag = 'Ca';
            app.element_20.Layout.Row = 4;
            app.element_20.Layout.Column = 2;
            app.element_20.ImageSource = '20.svg';

            % Create element_21
            app.element_21 = uiimage(app.GridLayout);
            app.element_21.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_21.Tag = 'Sc';
            app.element_21.Layout.Row = 4;
            app.element_21.Layout.Column = 3;
            app.element_21.ImageSource = '21.svg';

            % Create element_22
            app.element_22 = uiimage(app.GridLayout);
            app.element_22.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_22.Tag = 'Ti';
            app.element_22.Layout.Row = 4;
            app.element_22.Layout.Column = 4;
            app.element_22.ImageSource = '22.svg';

            % Create element_23
            app.element_23 = uiimage(app.GridLayout);
            app.element_23.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_23.Tag = 'V';
            app.element_23.Layout.Row = 4;
            app.element_23.Layout.Column = 5;
            app.element_23.ImageSource = '23.svg';

            % Create element_24
            app.element_24 = uiimage(app.GridLayout);
            app.element_24.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_24.Tag = 'Cr';
            app.element_24.Layout.Row = 4;
            app.element_24.Layout.Column = 6;
            app.element_24.ImageSource = '24.svg';

            % Create element_25
            app.element_25 = uiimage(app.GridLayout);
            app.element_25.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_25.Tag = 'Mn';
            app.element_25.Layout.Row = 4;
            app.element_25.Layout.Column = 7;
            app.element_25.ImageSource = '25.svg';

            % Create element_26
            app.element_26 = uiimage(app.GridLayout);
            app.element_26.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_26.Tag = 'Fe';
            app.element_26.Layout.Row = 4;
            app.element_26.Layout.Column = 8;
            app.element_26.ImageSource = '26.svg';

            % Create element_27
            app.element_27 = uiimage(app.GridLayout);
            app.element_27.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_27.Tag = 'Co';
            app.element_27.Layout.Row = 4;
            app.element_27.Layout.Column = 9;
            app.element_27.ImageSource = '27.svg';

            % Create element_28
            app.element_28 = uiimage(app.GridLayout);
            app.element_28.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_28.Tag = 'Ni';
            app.element_28.Layout.Row = 4;
            app.element_28.Layout.Column = 10;
            app.element_28.ImageSource = '28.svg';

            % Create element_29
            app.element_29 = uiimage(app.GridLayout);
            app.element_29.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_29.Tag = 'Cu';
            app.element_29.Layout.Row = 4;
            app.element_29.Layout.Column = 11;
            app.element_29.ImageSource = '29.svg';

            % Create element_30
            app.element_30 = uiimage(app.GridLayout);
            app.element_30.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_30.Tag = 'Zn';
            app.element_30.Layout.Row = 4;
            app.element_30.Layout.Column = 12;
            app.element_30.ImageSource = '30.svg';

            % Create element_31
            app.element_31 = uiimage(app.GridLayout);
            app.element_31.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_31.Tag = 'Ga';
            app.element_31.Layout.Row = 4;
            app.element_31.Layout.Column = 13;
            app.element_31.ImageSource = '31.svg';

            % Create element_32
            app.element_32 = uiimage(app.GridLayout);
            app.element_32.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_32.Tag = 'Ge';
            app.element_32.Layout.Row = 4;
            app.element_32.Layout.Column = 14;
            app.element_32.ImageSource = '32.svg';

            % Create element_33
            app.element_33 = uiimage(app.GridLayout);
            app.element_33.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_33.Tag = 'As';
            app.element_33.Layout.Row = 4;
            app.element_33.Layout.Column = 15;
            app.element_33.ImageSource = '33.svg';

            % Create element_34
            app.element_34 = uiimage(app.GridLayout);
            app.element_34.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_34.Tag = 'Se';
            app.element_34.Layout.Row = 4;
            app.element_34.Layout.Column = 16;
            app.element_34.ImageSource = '34.svg';

            % Create element_35
            app.element_35 = uiimage(app.GridLayout);
            app.element_35.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_35.Tag = 'Br';
            app.element_35.Layout.Row = 4;
            app.element_35.Layout.Column = 17;
            app.element_35.ImageSource = '35.svg';

            % Create element_36
            app.element_36 = uiimage(app.GridLayout);
            app.element_36.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_36.Tag = 'Kr';
            app.element_36.Layout.Row = 4;
            app.element_36.Layout.Column = 18;
            app.element_36.ImageSource = '36.svg';

            % Create element_37
            app.element_37 = uiimage(app.GridLayout);
            app.element_37.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_37.Tag = 'Rb';
            app.element_37.Layout.Row = 5;
            app.element_37.Layout.Column = 1;
            app.element_37.ImageSource = '37.svg';

            % Create element_38
            app.element_38 = uiimage(app.GridLayout);
            app.element_38.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_38.Tag = 'Sr';
            app.element_38.Layout.Row = 5;
            app.element_38.Layout.Column = 2;
            app.element_38.ImageSource = '38.svg';

            % Create element_39
            app.element_39 = uiimage(app.GridLayout);
            app.element_39.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_39.Tag = 'Y';
            app.element_39.Layout.Row = 5;
            app.element_39.Layout.Column = 3;
            app.element_39.ImageSource = '39.svg';

            % Create element_40
            app.element_40 = uiimage(app.GridLayout);
            app.element_40.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_40.Tag = 'Zr';
            app.element_40.Layout.Row = 5;
            app.element_40.Layout.Column = 4;
            app.element_40.ImageSource = '40.svg';

            % Create element_41
            app.element_41 = uiimage(app.GridLayout);
            app.element_41.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_41.Tag = 'Nb';
            app.element_41.Layout.Row = 5;
            app.element_41.Layout.Column = 5;
            app.element_41.ImageSource = '41.svg';

            % Create element_42
            app.element_42 = uiimage(app.GridLayout);
            app.element_42.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_42.Tag = 'Mo';
            app.element_42.Layout.Row = 5;
            app.element_42.Layout.Column = 6;
            app.element_42.ImageSource = '42.svg';

            % Create element_43
            app.element_43 = uiimage(app.GridLayout);
            app.element_43.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_43.Tag = 'Tc';
            app.element_43.Layout.Row = 5;
            app.element_43.Layout.Column = 7;
            app.element_43.ImageSource = '43.svg';

            % Create element_44
            app.element_44 = uiimage(app.GridLayout);
            app.element_44.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_44.Tag = 'Ru';
            app.element_44.Layout.Row = 5;
            app.element_44.Layout.Column = 8;
            app.element_44.ImageSource = '44.svg';

            % Create element_45
            app.element_45 = uiimage(app.GridLayout);
            app.element_45.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_45.Tag = 'Rh';
            app.element_45.Layout.Row = 5;
            app.element_45.Layout.Column = 9;
            app.element_45.ImageSource = '45.svg';

            % Create element_46
            app.element_46 = uiimage(app.GridLayout);
            app.element_46.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_46.Tag = 'Pd';
            app.element_46.Layout.Row = 5;
            app.element_46.Layout.Column = 10;
            app.element_46.ImageSource = '46.svg';

            % Create element_47
            app.element_47 = uiimage(app.GridLayout);
            app.element_47.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_47.Tag = 'Ag';
            app.element_47.Layout.Row = 5;
            app.element_47.Layout.Column = 11;
            app.element_47.ImageSource = '47.svg';

            % Create element_48
            app.element_48 = uiimage(app.GridLayout);
            app.element_48.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_48.Tag = 'CD';
            app.element_48.Layout.Row = 5;
            app.element_48.Layout.Column = 12;
            app.element_48.ImageSource = '48.svg';

            % Create element_49
            app.element_49 = uiimage(app.GridLayout);
            app.element_49.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_49.Tag = 'In';
            app.element_49.Layout.Row = 5;
            app.element_49.Layout.Column = 13;
            app.element_49.ImageSource = '49.svg';

            % Create element_50
            app.element_50 = uiimage(app.GridLayout);
            app.element_50.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_50.Tag = 'Sn';
            app.element_50.Layout.Row = 5;
            app.element_50.Layout.Column = 14;
            app.element_50.ImageSource = '50.svg';

            % Create element_51
            app.element_51 = uiimage(app.GridLayout);
            app.element_51.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_51.Tag = 'Sb';
            app.element_51.Layout.Row = 5;
            app.element_51.Layout.Column = 15;
            app.element_51.ImageSource = '51.svg';

            % Create element_52
            app.element_52 = uiimage(app.GridLayout);
            app.element_52.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_52.Tag = 'Te';
            app.element_52.Layout.Row = 5;
            app.element_52.Layout.Column = 16;
            app.element_52.ImageSource = '52.svg';

            % Create element_53
            app.element_53 = uiimage(app.GridLayout);
            app.element_53.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_53.Tag = 'I';
            app.element_53.Layout.Row = 5;
            app.element_53.Layout.Column = 17;
            app.element_53.ImageSource = '53.svg';

            % Create element_54
            app.element_54 = uiimage(app.GridLayout);
            app.element_54.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_54.Tag = 'Xe';
            app.element_54.Layout.Row = 5;
            app.element_54.Layout.Column = 18;
            app.element_54.ImageSource = '54.svg';

            % Create element_55
            app.element_55 = uiimage(app.GridLayout);
            app.element_55.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_55.Tag = 'Cs';
            app.element_55.Layout.Row = 6;
            app.element_55.Layout.Column = 1;
            app.element_55.ImageSource = '55.svg';

            % Create element_56
            app.element_56 = uiimage(app.GridLayout);
            app.element_56.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_56.Tag = 'Ba';
            app.element_56.Layout.Row = 6;
            app.element_56.Layout.Column = 2;
            app.element_56.ImageSource = '56.svg';

            % Create element_57_71
            app.element_57_71 = uiimage(app.GridLayout);
            app.element_57_71.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_57_71.Enable = 'off';
            app.element_57_71.Layout.Row = 6;
            app.element_57_71.Layout.Column = 3;
            app.element_57_71.ImageSource = '57_71.svg';

            % Create element_72
            app.element_72 = uiimage(app.GridLayout);
            app.element_72.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_72.Tag = 'Hf';
            app.element_72.Layout.Row = 6;
            app.element_72.Layout.Column = 4;
            app.element_72.ImageSource = '72.svg';

            % Create element_73
            app.element_73 = uiimage(app.GridLayout);
            app.element_73.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_73.Tag = 'Ta';
            app.element_73.Layout.Row = 6;
            app.element_73.Layout.Column = 5;
            app.element_73.ImageSource = '73.svg';

            % Create element_74
            app.element_74 = uiimage(app.GridLayout);
            app.element_74.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_74.Tag = 'W';
            app.element_74.Layout.Row = 6;
            app.element_74.Layout.Column = 6;
            app.element_74.ImageSource = '74.svg';

            % Create element_75
            app.element_75 = uiimage(app.GridLayout);
            app.element_75.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_75.Tag = 'Re';
            app.element_75.Layout.Row = 6;
            app.element_75.Layout.Column = 7;
            app.element_75.ImageSource = '75.svg';

            % Create element_76
            app.element_76 = uiimage(app.GridLayout);
            app.element_76.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_76.Tag = 'Os';
            app.element_76.Layout.Row = 6;
            app.element_76.Layout.Column = 8;
            app.element_76.ImageSource = '76.svg';

            % Create element_77
            app.element_77 = uiimage(app.GridLayout);
            app.element_77.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_77.Tag = 'Ir';
            app.element_77.Layout.Row = 6;
            app.element_77.Layout.Column = 9;
            app.element_77.ImageSource = '77.svg';

            % Create element_78
            app.element_78 = uiimage(app.GridLayout);
            app.element_78.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_78.Tag = 'Pt';
            app.element_78.Layout.Row = 6;
            app.element_78.Layout.Column = 10;
            app.element_78.ImageSource = '78.svg';

            % Create element_79
            app.element_79 = uiimage(app.GridLayout);
            app.element_79.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_79.Tag = 'Au';
            app.element_79.Layout.Row = 6;
            app.element_79.Layout.Column = 11;
            app.element_79.ImageSource = '79.svg';

            % Create element_80
            app.element_80 = uiimage(app.GridLayout);
            app.element_80.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_80.Tag = 'Hg';
            app.element_80.Layout.Row = 6;
            app.element_80.Layout.Column = 12;
            app.element_80.ImageSource = '80.svg';

            % Create element_81
            app.element_81 = uiimage(app.GridLayout);
            app.element_81.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_81.Tag = 'Tl';
            app.element_81.Layout.Row = 6;
            app.element_81.Layout.Column = 13;
            app.element_81.ImageSource = '81.svg';

            % Create element_82
            app.element_82 = uiimage(app.GridLayout);
            app.element_82.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_82.Tag = 'Pb';
            app.element_82.Layout.Row = 6;
            app.element_82.Layout.Column = 14;
            app.element_82.ImageSource = '82.svg';

            % Create element_83
            app.element_83 = uiimage(app.GridLayout);
            app.element_83.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_83.Tag = 'Bi';
            app.element_83.Layout.Row = 6;
            app.element_83.Layout.Column = 15;
            app.element_83.ImageSource = '83.svg';

            % Create element_84
            app.element_84 = uiimage(app.GridLayout);
            app.element_84.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_84.Tag = 'Po';
            app.element_84.Layout.Row = 6;
            app.element_84.Layout.Column = 16;
            app.element_84.ImageSource = '84.svg';

            % Create element_85
            app.element_85 = uiimage(app.GridLayout);
            app.element_85.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_85.Tag = 'At';
            app.element_85.Layout.Row = 6;
            app.element_85.Layout.Column = 17;
            app.element_85.ImageSource = '85.svg';

            % Create element_86
            app.element_86 = uiimage(app.GridLayout);
            app.element_86.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_86.Tag = 'Rn';
            app.element_86.Layout.Row = 6;
            app.element_86.Layout.Column = 18;
            app.element_86.ImageSource = '86.svg';

            % Create element_87
            app.element_87 = uiimage(app.GridLayout);
            app.element_87.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_87.Tag = 'Fr';
            app.element_87.Layout.Row = 7;
            app.element_87.Layout.Column = 1;
            app.element_87.ImageSource = '87.svg';

            % Create element_88
            app.element_88 = uiimage(app.GridLayout);
            app.element_88.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_88.Tag = 'Ra';
            app.element_88.Layout.Row = 7;
            app.element_88.Layout.Column = 2;
            app.element_88.ImageSource = '88.svg';

            % Create element_89
            app.element_89 = uiimage(app.GridLayout);
            app.element_89.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_89.Tag = 'Ac';
            app.element_89.Layout.Row = 7;
            app.element_89.Layout.Column = 3;
            app.element_89.ImageSource = '89.svg';

            % Create element_90
            app.element_90 = uiimage(app.GridLayout);
            app.element_90.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_90.Tag = 'Th';
            app.element_90.Layout.Row = 7;
            app.element_90.Layout.Column = 4;
            app.element_90.ImageSource = '90.svg';

            % Create element_91
            app.element_91 = uiimage(app.GridLayout);
            app.element_91.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_91.Tag = 'Pa';
            app.element_91.Layout.Row = 7;
            app.element_91.Layout.Column = 5;
            app.element_91.ImageSource = '91.svg';

            % Create element_92
            app.element_92 = uiimage(app.GridLayout);
            app.element_92.ImageClickedFcn = createCallbackFcn(app, @element_1ImageClicked, true);
            app.element_92.Tag = 'U';
            app.element_92.Layout.Row = 7;
            app.element_92.Layout.Column = 6;
            app.element_92.ImageSource = '92.svg';

            % Create Legend
            app.Legend = uiimage(app.GridLayout);
            app.Legend.Layout.Row = 1;
            app.Legend.Layout.Column = [13 17];
            app.Legend.ImageSource = 'legend_uielements.svg';

            % Create Database
            app.Database = uipanel(app.GridLayout);
            app.Database.BorderType = 'none';
            app.Database.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Database.Layout.Row = [1 3];
            app.Database.Layout.Column = [4 12];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.Database);
            app.GridLayout2.ColumnWidth = {'1.5x', 25, '1x', '1.5x', 25, '3.5x'};
            app.GridLayout2.RowHeight = {22, 26, 25, 23, 23, '3x'};
            app.GridLayout2.ColumnSpacing = 9.65625;
            app.GridLayout2.RowSpacing = 3.91294642857143;
            app.GridLayout2.Padding = [9.65625 3.91294642857143 9.65625 3.91294642857143];
            app.GridLayout2.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create text_list_species
            app.text_list_species = uilabel(app.GridLayout2);
            app.text_list_species.Layout.Row = 1;
            app.text_list_species.Layout.Column = [1 2];
            app.text_list_species.Text = 'Search Species';

            % Create Panel_5
            app.Panel_5 = uipanel(app.GridLayout2);
            app.Panel_5.BorderType = 'none';
            app.Panel_5.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_5.Layout.Row = [1 6];
            app.Panel_5.Layout.Column = 6;

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.Panel_5);
            app.GridLayout3.ColumnWidth = {70, 71, '2.66x'};
            app.GridLayout3.RowHeight = {19, 19, 19, 46, 19, 19, 19};
            app.GridLayout3.RowSpacing = 8.0703125;
            app.GridLayout3.Padding = [10 8.0703125 0 8.0703125];
            app.GridLayout3.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create text_species_label
            app.text_species_label = uilabel(app.GridLayout3);
            app.text_species_label.HorizontalAlignment = 'center';
            app.text_species_label.Layout.Row = 1;
            app.text_species_label.Layout.Column = 1;
            app.text_species_label.Text = 'Species';

            % Create text_species
            app.text_species = uieditfield(app.GridLayout3, 'text');
            app.text_species.Editable = 'off';
            app.text_species.HorizontalAlignment = 'center';
            app.text_species.Layout.Row = 1;
            app.text_species.Layout.Column = [2 3];

            % Create text_codename_label
            app.text_codename_label = uilabel(app.GridLayout3);
            app.text_codename_label.HorizontalAlignment = 'center';
            app.text_codename_label.Layout.Row = 2;
            app.text_codename_label.Layout.Column = 1;
            app.text_codename_label.Text = 'Codename';

            % Create text_codename
            app.text_codename = uieditfield(app.GridLayout3, 'text');
            app.text_codename.Editable = 'off';
            app.text_codename.HorizontalAlignment = 'center';
            app.text_codename.Layout.Row = 2;
            app.text_codename.Layout.Column = [2 3];

            % Create text_phase_label
            app.text_phase_label = uilabel(app.GridLayout3);
            app.text_phase_label.HorizontalAlignment = 'center';
            app.text_phase_label.Layout.Row = 3;
            app.text_phase_label.Layout.Column = 1;
            app.text_phase_label.Text = 'Phase';

            % Create text_phase
            app.text_phase = uieditfield(app.GridLayout3, 'text');
            app.text_phase.Editable = 'off';
            app.text_phase.HorizontalAlignment = 'center';
            app.text_phase.Layout.Row = 3;
            app.text_phase.Layout.Column = [2 3];

            % Create text_W_label
            app.text_W_label = uilabel(app.GridLayout3);
            app.text_W_label.HorizontalAlignment = 'center';
            app.text_W_label.Layout.Row = 5;
            app.text_W_label.Layout.Column = [1 2];
            app.text_W_label.Text = 'Molecular Weight [g/mol]';

            % Create text_W
            app.text_W = uieditfield(app.GridLayout3, 'numeric');
            app.text_W.ValueDisplayFormat = '%.6g';
            app.text_W.Editable = 'off';
            app.text_W.HorizontalAlignment = 'center';
            app.text_W.Layout.Row = 5;
            app.text_W.Layout.Column = 3;

            % Create text_hf_label
            app.text_hf_label = uilabel(app.GridLayout3);
            app.text_hf_label.HorizontalAlignment = 'center';
            app.text_hf_label.Layout.Row = 6;
            app.text_hf_label.Layout.Column = [1 2];
            app.text_hf_label.Text = 'Enthalpy formation [kJ]';

            % Create text_hf
            app.text_hf = uieditfield(app.GridLayout3, 'numeric');
            app.text_hf.ValueDisplayFormat = '%.6g';
            app.text_hf.Editable = 'off';
            app.text_hf.HorizontalAlignment = 'center';
            app.text_hf.Layout.Row = 6;
            app.text_hf.Layout.Column = 3;

            % Create text_ef_label
            app.text_ef_label = uilabel(app.GridLayout3);
            app.text_ef_label.HorizontalAlignment = 'center';
            app.text_ef_label.Layout.Row = 7;
            app.text_ef_label.Layout.Column = [1 2];
            app.text_ef_label.Text = 'Int. energy formation [kJ]';

            % Create text_ef
            app.text_ef = uieditfield(app.GridLayout3, 'numeric');
            app.text_ef.ValueDisplayFormat = '%.6g';
            app.text_ef.Editable = 'off';
            app.text_ef.HorizontalAlignment = 'center';
            app.text_ef.Layout.Row = 7;
            app.text_ef.Layout.Column = 3;

            % Create CommentsTextAreaLabel
            app.CommentsTextAreaLabel = uilabel(app.GridLayout3);
            app.CommentsTextAreaLabel.HorizontalAlignment = 'center';
            app.CommentsTextAreaLabel.Layout.Row = 4;
            app.CommentsTextAreaLabel.Layout.Column = 1;
            app.CommentsTextAreaLabel.Text = 'Comments';

            % Create text_comments
            app.text_comments = uitextarea(app.GridLayout3);
            app.text_comments.Layout.Row = 4;
            app.text_comments.Layout.Column = [2 3];

            % Create edit_seeker
            app.edit_seeker = uieditfield(app.GridLayout2, 'text');
            app.edit_seeker.ValueChangingFcn = createCallbackFcn(app, @edit_seekerValueChanging, true);
            app.edit_seeker.Layout.Row = 2;
            app.edit_seeker.Layout.Column = [1 2];

            % Create text_LS_DB
            app.text_LS_DB = uilabel(app.GridLayout2);
            app.text_LS_DB.Layout.Row = 3;
            app.text_LS_DB.Layout.Column = [1 2];
            app.text_LS_DB.Text = 'List of species DB';

            % Create text_LS
            app.text_LS = uilabel(app.GridLayout2);
            app.text_LS.Layout.Row = 3;
            app.text_LS.Layout.Column = [4 5];
            app.text_LS.Text = 'List of species (export)';

            % Create Copy
            app.Copy = uibutton(app.GridLayout2, 'push');
            app.Copy.ButtonPushedFcn = createCallbackFcn(app, @CopyButtonPushed, true);
            app.Copy.Icon = 'icon_copy.svg';
            app.Copy.IconAlignment = 'center';
            app.Copy.Tooltip = {''};
            app.Copy.Layout.Row = 3;
            app.Copy.Layout.Column = 5;
            app.Copy.Text = '';

            % Create listbox_LS_DB
            app.listbox_LS_DB = uilistbox(app.GridLayout2);
            app.listbox_LS_DB.Items = {};
            app.listbox_LS_DB.Multiselect = 'on';
            app.listbox_LS_DB.ValueChangedFcn = createCallbackFcn(app, @listbox_LSValueChanged, true);
            app.listbox_LS_DB.Layout.Row = [4 6];
            app.listbox_LS_DB.Layout.Column = [1 2];
            app.listbox_LS_DB.ClickedFcn = createCallbackFcn(app, @listbox_LS_DBClicked, true);
            app.listbox_LS_DB.Value = {};

            % Create AddButton1
            app.AddButton1 = uibutton(app.GridLayout2, 'push');
            app.AddButton1.ButtonPushedFcn = createCallbackFcn(app, @AddButton1Pushed, true);
            app.AddButton1.Layout.Row = 4;
            app.AddButton1.Layout.Column = 3;
            app.AddButton1.Text = 'Add >>';

            % Create listbox_LS
            app.listbox_LS = uilistbox(app.GridLayout2);
            app.listbox_LS.Items = {};
            app.listbox_LS.Multiselect = 'on';
            app.listbox_LS.ValueChangedFcn = createCallbackFcn(app, @listbox_LSValueChanged, true);
            app.listbox_LS.Layout.Row = [4 6];
            app.listbox_LS.Layout.Column = [4 5];
            app.listbox_LS.ClickedFcn = createCallbackFcn(app, @listbox_LSClicked, true);
            app.listbox_LS.Value = {};

            % Create RemoveButton1
            app.RemoveButton1 = uibutton(app.GridLayout2, 'push');
            app.RemoveButton1.ButtonPushedFcn = createCallbackFcn(app, @RemoveButton1Pushed, true);
            app.RemoveButton1.Layout.Row = 5;
            app.RemoveButton1.Layout.Column = 3;
            app.RemoveButton1.Text = 'Remove <<';

            % Create plot
            app.plot = uibutton(app.GridLayout2, 'push');
            app.plot.ButtonPushedFcn = createCallbackFcn(app, @plotButtonPushed, true);
            app.plot.Layout.Row = 1;
            app.plot.Layout.Column = 3;
            app.plot.Text = 'Plot';

            % Create property
            app.property = uidropdown(app.GridLayout2);
            app.property.Items = {'Specific heat capacity at constant pressure', 'Specific heat capacity at constant volume', 'Adiabatic index', 'Internal energy', 'Enthalpy', 'Entropy', 'Gibbs energy'};
            app.property.ItemsData = {'cp', 'cv', 'gamma', 'e', 'h', 's', 'g'};
            app.property.Layout.Row = 2;
            app.property.Layout.Column = [4 5];
            app.property.Value = 'cp';

            % Create temperature_range
            app.temperature_range = uieditfield(app.GridLayout2, 'text');
            app.temperature_range.HorizontalAlignment = 'center';
            app.temperature_range.Tooltip = {'Select temperature range [initial:step:final] values'};
            app.temperature_range.Layout.Row = 1;
            app.temperature_range.Layout.Column = 4;
            app.temperature_range.Value = '[300:10:20000]';

            % Create type
            app.type = uidropdown(app.GridLayout2);
            app.type.Items = {'CT', 'NASA'};
            app.type.Tooltip = {'Select polynomial fits'};
            app.type.Layout.Row = 2;
            app.type.Layout.Column = 3;
            app.type.Value = 'CT';

            % Create log
            app.log = uicheckbox(app.GridLayout2);
            app.log.Tooltip = {'Set logarithmic scale'};
            app.log.Text = '';
            app.log.Layout.Row = 1;
            app.log.Layout.Column = 5;

            % Create Panel_Export
            app.Panel_Export = uipanel(app.GridLayout);
            app.Panel_Export.BorderType = 'none';
            app.Panel_Export.TitlePosition = 'centertop';
            app.Panel_Export.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_Export.Layout.Row = 7;
            app.Panel_Export.Layout.Column = 17;

            % Create Export
            app.Export = uibutton(app.Panel_Export, 'push');
            app.Export.ButtonPushedFcn = createCallbackFcn(app, @ExportButtonPushed, true);
            app.Export.IconAlignment = 'center';
            app.Export.Position = [1 2 69 25];
            app.Export.Text = 'Export';

            % Create Panel_Close
            app.Panel_Close = uipanel(app.GridLayout);
            app.Panel_Close.BorderType = 'none';
            app.Panel_Close.TitlePosition = 'centertop';
            app.Panel_Close.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_Close.Layout.Row = 7;
            app.Panel_Close.Layout.Column = 18;

            % Create Close
            app.Close = uibutton(app.Panel_Close, 'push');
            app.Close.ButtonPushedFcn = createCallbackFcn(app, @UIElementsCloseRequest, true);
            app.Close.IconAlignment = 'center';
            app.Close.Position = [1 2 69 25];
            app.Close.Text = 'Close';

            % Create Panel_6
            app.Panel_6 = uipanel(app.GridLayout);
            app.Panel_6.BorderType = 'none';
            app.Panel_6.BackgroundColor = [0.9098 0.9098 0.8902];
            app.Panel_6.Layout.Row = [2 3];
            app.Panel_6.Layout.Column = 3;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.Panel_6);
            app.GridLayout4.ColumnWidth = {84};
            app.GridLayout4.RowHeight = {22, 22, 22, '1x'};
            app.GridLayout4.RowSpacing = 1;
            app.GridLayout4.Padding = [0 10 0 10];
            app.GridLayout4.BackgroundColor = [0.9098 0.9098 0.8902];

            % Create ionized
            app.ionized = uicheckbox(app.GridLayout4);
            app.ionized.ValueChangedFcn = createCallbackFcn(app, @ionizedValueChanged, true);
            app.ionized.Tooltip = {'Consider ionized species?'};
            app.ionized.Text = 'Ionized';
            app.ionized.Layout.Row = 1;
            app.ionized.Layout.Column = 1;

            % Create burcat
            app.burcat = uicheckbox(app.GridLayout4);
            app.burcat.ValueChangedFcn = createCallbackFcn(app, @burcatValueChanged, true);
            app.burcat.Tooltip = {'Consider species from Third Millennium database (Burcat)?'};
            app.burcat.Text = 'Burcat';
            app.burcat.Layout.Row = 3;
            app.burcat.Layout.Column = 1;

            % Create condensed
            app.condensed = uicheckbox(app.GridLayout4);
            app.condensed.ValueChangedFcn = createCallbackFcn(app, @condensedValueChanged, true);
            app.condensed.Tooltip = {'Consider condensed species?'};
            app.condensed.Text = 'Conden.';
            app.condensed.Layout.Row = 2;
            app.condensed.Layout.Column = 1;

            % Create ContextMenu
            app.ContextMenu = uicontextmenu(app.UIElements);

            % Create enableremoveMenu
            app.enableremoveMenu = uimenu(app.ContextMenu);
            app.enableremoveMenu.MenuSelectedFcn = createCallbackFcn(app, @enableremoveMenuSelected, true);
            app.enableremoveMenu.Visible = 'off';
            app.enableremoveMenu.Text = 'enable/remove';
            
            % Assign app.ContextMenu
            app.element_1.ContextMenu = app.ContextMenu;
            app.element_1_2.ContextMenu = app.ContextMenu;
            app.element_1_3.ContextMenu = app.ContextMenu;
            app.element_2.ContextMenu = app.ContextMenu;
            app.element_3.ContextMenu = app.ContextMenu;
            app.element_4.ContextMenu = app.ContextMenu;
            app.element_5.ContextMenu = app.ContextMenu;
            app.element_6.ContextMenu = app.ContextMenu;
            app.element_7.ContextMenu = app.ContextMenu;
            app.element_8.ContextMenu = app.ContextMenu;
            app.element_9.ContextMenu = app.ContextMenu;
            app.element_10.ContextMenu = app.ContextMenu;
            app.element_11.ContextMenu = app.ContextMenu;
            app.element_12.ContextMenu = app.ContextMenu;
            app.element_13.ContextMenu = app.ContextMenu;
            app.element_14.ContextMenu = app.ContextMenu;
            app.element_15.ContextMenu = app.ContextMenu;
            app.element_16.ContextMenu = app.ContextMenu;
            app.element_17.ContextMenu = app.ContextMenu;
            app.element_18.ContextMenu = app.ContextMenu;
            app.element_19.ContextMenu = app.ContextMenu;
            app.element_20.ContextMenu = app.ContextMenu;
            app.element_21.ContextMenu = app.ContextMenu;
            app.element_22.ContextMenu = app.ContextMenu;
            app.element_23.ContextMenu = app.ContextMenu;
            app.element_24.ContextMenu = app.ContextMenu;
            app.element_25.ContextMenu = app.ContextMenu;
            app.element_26.ContextMenu = app.ContextMenu;
            app.element_27.ContextMenu = app.ContextMenu;
            app.element_28.ContextMenu = app.ContextMenu;
            app.element_29.ContextMenu = app.ContextMenu;
            app.element_30.ContextMenu = app.ContextMenu;
            app.element_31.ContextMenu = app.ContextMenu;
            app.element_32.ContextMenu = app.ContextMenu;
            app.element_33.ContextMenu = app.ContextMenu;
            app.element_34.ContextMenu = app.ContextMenu;
            app.element_35.ContextMenu = app.ContextMenu;
            app.element_36.ContextMenu = app.ContextMenu;
            app.element_37.ContextMenu = app.ContextMenu;
            app.element_38.ContextMenu = app.ContextMenu;
            app.element_39.ContextMenu = app.ContextMenu;
            app.element_40.ContextMenu = app.ContextMenu;
            app.element_41.ContextMenu = app.ContextMenu;
            app.element_42.ContextMenu = app.ContextMenu;
            app.element_43.ContextMenu = app.ContextMenu;
            app.element_44.ContextMenu = app.ContextMenu;
            app.element_45.ContextMenu = app.ContextMenu;
            app.element_46.ContextMenu = app.ContextMenu;
            app.element_47.ContextMenu = app.ContextMenu;
            app.element_48.ContextMenu = app.ContextMenu;
            app.element_49.ContextMenu = app.ContextMenu;
            app.element_50.ContextMenu = app.ContextMenu;
            app.element_51.ContextMenu = app.ContextMenu;
            app.element_52.ContextMenu = app.ContextMenu;
            app.element_53.ContextMenu = app.ContextMenu;
            app.element_54.ContextMenu = app.ContextMenu;
            app.element_55.ContextMenu = app.ContextMenu;
            app.element_56.ContextMenu = app.ContextMenu;
            app.element_72.ContextMenu = app.ContextMenu;
            app.element_73.ContextMenu = app.ContextMenu;
            app.element_74.ContextMenu = app.ContextMenu;
            app.element_75.ContextMenu = app.ContextMenu;
            app.element_76.ContextMenu = app.ContextMenu;
            app.element_77.ContextMenu = app.ContextMenu;
            app.element_78.ContextMenu = app.ContextMenu;
            app.element_79.ContextMenu = app.ContextMenu;
            app.element_80.ContextMenu = app.ContextMenu;
            app.element_81.ContextMenu = app.ContextMenu;
            app.element_82.ContextMenu = app.ContextMenu;
            app.element_83.ContextMenu = app.ContextMenu;
            app.element_84.ContextMenu = app.ContextMenu;
            app.element_85.ContextMenu = app.ContextMenu;
            app.element_86.ContextMenu = app.ContextMenu;
            app.element_87.ContextMenu = app.ContextMenu;
            app.element_88.ContextMenu = app.ContextMenu;
            app.element_89.ContextMenu = app.ContextMenu;
            app.element_90.ContextMenu = app.ContextMenu;
            app.element_91.ContextMenu = app.ContextMenu;
            app.element_92.ContextMenu = app.ContextMenu;

            % Show the figure after all components are created
            app.UIElements.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = uielements(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.UIElements)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.UIElements)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIElements)
        end
    end
end
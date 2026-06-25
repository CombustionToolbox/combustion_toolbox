classdef AppController < handle
    % Coordinates GUI callbacks, state, panels, solving, and results.
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %     session (AppSession): Long-lived CT services used by the GUI
    %     state (AppState): Current GUI state without App Designer handles
    %     input (AppInput): Last typed input snapshot read from the GUI
    %     solver (AppSolver): GUI-level solver orchestration object
    %     results (AppResults): GUI-side storage for computed studies
    %     problemCatalog (AppProblemCatalog): Problem definitions and metadata
    %     problemPanel (AppProblemPanel): Problem-selection UI updater
    %     problemInputsPanel (AppProblemInputsPanel): Coupled input updater
    %     speciesPanel (AppSpeciesPanel): Reactants/products UI updater
    %     resultsPanel (AppResultsPanel): Results tree, tables, and status updater
    %     statusPanel (AppStatusPanel): Status indicator updater
    %     consolePanel (AppConsolePanel): Command-window updater
    %     currentSolution (struct): Last successful solver solution
    %
    % Examples:
    %     * controller = combustiontoolbox.gui.AppController(app);
    %     * input = controller.readInput();
    %     * definition = controller.currentProblem();

    properties (Access = private)
        app
        session
        state
        input
        solver
        results
        problemCatalog
        problemPanel
        problemInputsPanel
        speciesPanel
        resultsPanel
        statusPanel
        consolePanel
        currentSolution = []
    end

    methods
        function obj = AppController(app, varargin)
            % AppController constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Optional Args:
            %     session (AppSession): Existing session object. If omitted,
            %         the constructor references services already owned by app
            %
            % Returns:
            %     obj (AppController): Initialized controller object

            if nargin == 0
                app = [];
            end

            obj.app = app;

            if nargin > 1 && isa(varargin{1}, 'combustiontoolbox.gui.AppSession')
                obj.session = varargin{1};
            elseif nargin > 0 && ~isempty(app)
                obj.session = combustiontoolbox.gui.AppSession.fromApp(app);
            else
                obj.session = combustiontoolbox.gui.AppSession('initialize', false);
            end

            obj.input = combustiontoolbox.gui.AppInput();
            obj.state = combustiontoolbox.gui.AppState();
            obj.results = combustiontoolbox.gui.AppResults();
            obj.problemCatalog = combustiontoolbox.gui.AppProblemCatalog();
            obj.problemInputsPanel = combustiontoolbox.gui.AppProblemInputsPanel(app);
            obj.problemPanel = combustiontoolbox.gui.AppProblemPanel(app, obj.problemCatalog, obj.problemInputsPanel);
            obj.speciesPanel = combustiontoolbox.gui.AppSpeciesPanel(app);
            obj.resultsPanel = combustiontoolbox.gui.AppResultsPanel(app);

            if ~isempty(app)
                obj.statusPanel = app.statusPanel;
                obj.consolePanel = app.consolePanel;
            end

            obj.solver = combustiontoolbox.gui.AppSolver(obj.session, obj.problemCatalog);
        end

        function input = readInput(obj)
            % Read the current GUI state into an AppInput object
            %
            % Returns:
            %     input (AppInput): Typed snapshot of the current GUI inputs
            input = combustiontoolbox.gui.AppInput.fromApp(obj.app);
            obj.input = input;
            obj.state.updateFromInput(input);
        end

        function definition = currentProblem(obj)
            % Get the catalog definition for the currently selected problem
            %
            % Returns:
            %     definition (struct): Problem metadata from AppProblemCatalog
            input = readInput(obj);
            definition = obj.problemCatalog.get(input.problemType);
        end

        function [items, itemsData] = problemDropdownData(obj)
            % Return labels and values for the problem dropdown
            %
            % Returns:
            %     items (cell): Human-readable dropdown labels
            %     itemsData (cell): Problem identifiers used as dropdown values
            [items, itemsData] = obj.problemCatalog.dropdownData();
        end

        function onConsoleValueChanged(obj, event) %#ok<INUSD>
            % Run a command from the GUI command window
            
            if isempty(obj.app)
                return
            end

            commands = obj.consoleCommands();

            if isempty(commands)
                return
            end

            command = commands{1};
            FLAG_ERROR = false;

            try
                output = obj.runConsoleCommand(command, commands);
            catch exception
                output = obj.formatConsoleError(exception);
                obj.setStatus('error');
                FLAG_ERROR = true;
            end

            if ~isempty(output) || strcmp(command, 'clc')
                obj.consolePanel.setOutput(output);
            end

            obj.consolePanel.clearInput();

            if ~FLAG_ERROR
                obj.setStatus('idle');
            end
        end

        function onSnapshotSelected(obj)
            % Export the current app figure as an image or PDF
            filter = {'*.pdf'; '*.jpg'; '*.png'; '*.tif'};
            [filename, filepath] = uiputfile(filter);

            if ischar(filename)
                exportapp(obj.app.UIFigure, [filepath filename]);
            end
        end

        function definition = problemDefinition(obj, varargin)
            % Get catalog metadata for a selected or explicit problem
            %
            % Optional Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     definition (struct): Problem metadata from AppProblemCatalog
            if nargin > 1 && ~isempty(varargin{1})
                definition = obj.problemCatalog.get(varargin{1});
                return
            end

            definition = obj.currentProblem();
        end

        function value = problemDefaults(obj, varargin)
            % Get default GUI input values for a selected or explicit problem
            %
            % Optional Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Default values keyed by GUI component name
            definition = obj.problemDefinition(varargin{:});
            value = definition.defaultInputs;
        end

        function value = problemVisibility(obj, varargin)
            % Get visibility metadata for a selected or explicit problem
            %
            % Optional Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Logical visibility flags
            definition = obj.problemDefinition(varargin{:});
            value = definition.visibility;
        end

        function value = problemLabels(obj, varargin)
            % Get GUI labels for a selected or explicit problem
            %
            % Optional Args:
            %     problemType (char): Problem type identifier
            %
            % Returns:
            %     value (struct): Label text keyed by GUI component name
            definition = obj.problemDefinition(varargin{:});
            value = definition.labels;
        end

        function definition = onProblemTypeChanged(obj)
            % Apply selected problem metadata and refresh typed input
            %
            % Returns:
            %     definition (struct): Applied problem definition
            definition = obj.problemPanel.applySelectedProblem();

            if ~isempty(obj.app)
                obj.readInput();
            end
        end

        function solution = solveCurrent(obj)
            % Solve the current GUI study and store normalized results
            %
            % Returns:
            %     solution (struct): Solver output and GUI result metadata
            input = obj.readInput();
            currentMixture = [];

            if ~isempty(obj.app) && isprop(obj.app, 'mixture')
                currentMixture = obj.app.mixture;
            elseif ~isempty(obj.session)
                currentMixture = obj.session.mixture;
            end

            if ~isempty(obj.session)
                obj.session.mixture = currentMixture;
            end

            solution = obj.solver.solve(input, currentMixture);
            obj.currentSolution = solution;
            obj.results.set(solution.results, solution.results);
            obj.results.setLayout(solution.resultLayout);

            if ~isempty(obj.resultsPanel)
                obj.resultsPanel.applySolution(solution);
            end
        end

        function solution = onCalculatePushed(obj)
            % Solve the current problem from the calculate button
            solution = [];
            obj.setWorkingStatus();

            try
                solution = obj.solveCurrent();
                obj.runAutoReport(solution);
                obj.setFinishedStatus();
            catch exception
                obj.setErrorStatus(exception);
            end
        end

        function value = hasResults(obj)
            % Check whether solved results are stored
            %
            % Returns:
            %     value (logical): True if solved results are available
            value = ~isempty(obj.results) && obj.results.hasResults();
        end

        function solution = currentSolutionData(obj)
            % Return the last successful solver solution
            %
            % Returns:
            %     solution (struct): Last successful solver solution
            solution = obj.currentSolution;
        end

        function loadSolution(obj, solution)
            % Load a solver solution into controller state
            %
            % Args:
            %     solution (struct): Solver solution to store
            obj.currentSolution = solution;
        end

        function syncSessionFromApp(obj)
            % Synchronize session services from the App Designer shell
            if isempty(obj.app) || isempty(obj.session)
                return
            end

            if isprop(obj.app, 'chemicalSystem')
                obj.session.chemicalSystem = obj.app.chemicalSystem;
            end

            if isprop(obj.app, 'mixture')
                obj.session.mixture = obj.app.mixture;
            end
        end

        function onCustomFigurePlotPushed(obj)
            % Plot the selected custom figure variables
            obj.resultsPanel.plotCustomFigure();
        end

        function onExportSelected(obj, format)
            % Export the current semantic result states
            obj.resultsPanel.exportResults(format);
        end

        function onClearSelected(obj)
            % Clear current results and visible result controls
            obj.clearCurrentStudy();
        end

        function onNewSelected(obj)
            % Start a new app study from the current visible inputs
            obj.clearCurrentStudy();
        end

        function [filepath, script] = onExportProblemScriptSelected(obj)
            % Export the current problem setup to a MATLAB script
            filepath = [];
            script = '';

            try
                [filepath, script] = obj.exportCurrentProblemScript();
                obj.setScriptExportStatus(filepath);
            catch exception
                obj.setErrorStatus(exception);
            end
        end

        function [filepath, script] = exportCurrentProblemScript(obj, varargin)
            % Export the current problem setup to a MATLAB script file
            %
            % Optional Args:
            %     filepath (char): Output file path
            %
            % Returns:
            %     filepath (char): Written file path, or [] if cancelled
            %     script (char): Generated MATLAB script text
            setup = obj.currentProblemSetup();
            script = obj.solver.exportScript(setup);

            if nargin > 1
                filepath = varargin{1};
            else
                filepath = obj.selectProblemScriptFile(setup);
            end

            if isempty(filepath)
                return
            end

            filepath = obj.ensureMatlabScriptExtension(filepath);
            obj.writeTextFile(filepath, script);
        end

        function onResultsTreeSelectionChanged(obj)
            % Update result fields from the selected results tree node
            obj.resultsPanel.onTreeSelectionChanged();
        end

        function onReactantsChanged(obj, event)
            % Update mixture state after reactant selection
            %
            % Args:
            %     event (object): App Designer callback event
            obj.speciesPanel.onReactantsChanged(event);
            obj.problemInputsPanel.updateMachOrVelocity('Mach');
            obj.readInput();
        end

        function onProductsChanged(obj)
            % Update product species and dependent display lists
            obj.speciesPanel.onProductsChanged();
            obj.readInput();
        end

        function onReactantsTableEdited(obj, event)
            % Rebuild mixture state after reactants table edits
            %
            % Args:
            %     event (object): App Designer callback event
            obj.speciesPanel.onReactantsTableEdited(event);
            obj.problemInputsPanel.updateMachOrVelocity('Mach');
            obj.readInput();
        end

        function onEquivalenceRatioChanged(obj, event)
            % Recompute reactant composition for a new equivalence ratio
            %
            % Args:
            %     event (object): App Designer callback event
            obj.speciesPanel.onEquivalenceRatioChanged(event);
            obj.readInput();
        end

        function onReactantsTemperatureChanged(obj)
            % Apply reactant temperature and refresh shock inputs
            obj.speciesPanel.onReactantsTemperatureChanged();
            obj.problemInputsPanel.updateMachOrVelocity('Mach');
            obj.readInput();
        end

        function onFrozenChemistryChanged(obj)
            % Refresh caloric model and product species
            obj.speciesPanel.onFrozenChemistryChanged();
            obj.readInput();
        end

        function onIonizedSpeciesChanged(obj)
            % Refresh product species after ionization changes
            obj.speciesPanel.onIonizedSpeciesChanged();
            obj.readInput();
        end

        function onRocketModelChanged(obj)
            % Apply rocket model controls
            obj.problemPanel.onRocketModelChanged();
            obj.readInput();
        end

        function onProductSpeciesAdded(obj)
            % Add selected species to the product species list
            obj.speciesPanel.onProductSpeciesAdded();
            obj.readInput();
        end

        function onProductSpeciesRemoved(obj)
            % Remove selected species from the product species list
            obj.speciesPanel.onProductSpeciesRemoved();
            obj.readInput();
        end

        function onDatabaseSpeciesSearchChanging(obj, event)
            % Filter the database species list
            obj.speciesPanel.onDatabaseSpeciesSearchChanging(event);
        end

        function onDisplaySpeciesSearchChanging(obj, event)
            % Filter the display species lists
            obj.speciesPanel.onDisplaySpeciesSearchChanging(event);
        end

        function onDisplaySpeciesAdded(obj)
            % Add selected species to the display species list
            obj.speciesPanel.onDisplaySpeciesAdded();
            obj.readInput();
        end

        function onDisplaySpeciesRemoved(obj)
            % Remove selected species from the display species list
            obj.speciesPanel.onDisplaySpeciesRemoved();
            obj.readInput();
        end

        function onProductPressureChanged(obj, event)
            % Synchronize product pressure edits
            %
            % Args:
            %     event (object): App Designer callback event
            obj.problemInputsPanel.onProductPressureChanged(event);
            obj.readInput();
        end

        function onReactantPressureChanged(obj, event)
            % Synchronize reactant pressure edits and shock inputs
            %
            % Args:
            %     event (object): App Designer callback event
            obj.problemInputsPanel.onReactantPressureChanged(event);
            obj.readInput();
        end

        function onReactantMachChanged(obj)
            % Refresh shock velocity after Mach number edits
            obj.problemInputsPanel.onReactantMachChanged();
            obj.readInput();
        end

        function onReactantVelocityChanged(obj)
            % Refresh shock Mach number after velocity edits
            obj.problemInputsPanel.onReactantVelocityChanged();
            obj.readInput();
        end

        function onRocketReactantAreaChanging(obj)
            % Keep rocket area ratio inputs mutually exclusive
            obj.problemInputsPanel.onRocketReactantAreaChanging();
        end

        function onRocketProductAreaChanging(obj)
            % Keep rocket area ratio inputs mutually exclusive
            obj.problemInputsPanel.onRocketProductAreaChanging();
        end

        function onDetonationReactantAngleChanging(obj)
            % Keep detonation angle inputs mutually exclusive
            obj.problemInputsPanel.onDetonationReactantAngleChanging();
        end

        function onDetonationProductAngleChanging(obj)
            % Keep detonation angle inputs mutually exclusive
            obj.problemInputsPanel.onDetonationProductAngleChanging();
        end

        function onShockReactantAngleChanging(obj)
            % Keep shock angle inputs mutually exclusive
            obj.problemInputsPanel.onShockReactantAngleChanging();
        end

        function onShockProductAngleChanging(obj)
            % Keep shock angle inputs mutually exclusive
            obj.problemInputsPanel.onShockProductAngleChanging();
        end

    end

    methods (Access = private)
        function clearCurrentStudy(obj)
            obj.currentSolution = [];
            obj.input = combustiontoolbox.gui.AppInput();
            obj.state = combustiontoolbox.gui.AppState();

            if ~isempty(obj.results)
                obj.results.clear();
            end

            if ~isempty(obj.resultsPanel)
                obj.resultsPanel.clear();
            end

            obj.clearAppStatus();
        end

        function clearAppStatus(obj)
            if isempty(obj.app)
                return
            end

            obj.consolePanel.clearOutput();

            obj.setStatus('idle');
        end

        function setWorkingStatus(obj)
            if isempty(obj.app)
                return
            end

            obj.setStatus('working');
            obj.consolePanel.setOutput('Solving problem...');
        end

        function runAutoReport(obj, solution)
            if isempty(solution) || ~isfield(solution, 'solver') || ~isfield(solution, 'reportOutputs')
                return
            end

            if isempty(obj.input) || ~strcmpi(char(obj.input.reportType), 'Auto')
                return
            end

            if ismethod(solution.solver, 'report')
                solution.solver.report(solution.reportOutputs{:});
            end
        end

        function commands = consoleCommands(obj)
            if ~isempty(obj.consolePanel)
                commands = obj.consolePanel.inputCommands();
                return
            end

            commands = {};
        end

        function output = runConsoleCommand(obj, command, commands)
            command = char(command);

            switch lower(command)
                case 'about'
                    uiabout();
                    output = 'Running uiabout...';
                case 'clear'
                    obj.consolePanel.clearOutput();
                    output = ' ';
                case {'docs', 'documentation'}
                    combustiontoolbox.utils.SystemUtils.openWebsite( ...
                        combustiontoolbox.utils.SystemUtils.url.docs);
                    output = 'Redirecting to CT documentation using the default browser...';
                case 'license'
                    uilicense();
                    output = 'Running uilicense...';
                case 'update'
                    [~, output] = combustiontoolbox.utils.checkUpdate(obj.app.UIFigure);
                case 'validations'
                    uivalidations();
                    output = 'Running uivalidations...';
                case 'version'
                    output = sprintf('Version: %s\nDate: %s', ...
                        combustiontoolbox.common.Constants.release, ...
                        combustiontoolbox.common.Constants.date);
                case {'save', 'export', 'export(''xls'')', 'export(''csv'')'}
                    obj.onExportSelected('.xls');
                    output = 'Exporting results...';
                case 'export(''mat'')'
                    obj.onExportSelected('.mat');
                    output = 'Exporting results...';
                case {'export(''m'')', 'export(''script'')', 'script'}
                    obj.onExportProblemScriptSelected();
                    output = 'Exporting problem script...';
                case {'settings', 'configuration', 'uipreferences'}
                    uipreferences(obj.app);
                    output = 'Opening uipreferences...';
                case {'web', 'website', 'webpage'}
                    combustiontoolbox.utils.SystemUtils.websiteCT;
                    output = 'Redirecting to CT website using the default browser...';
                otherwise
                    output = obj.runConsoleEval(command, commands);
            end
        end

        function output = runConsoleEval(obj, command, commands)
            obj.consolePanel.setOutput(sprintf('Running: %s...', command));
            obj.setStatus('working');
            pause(0.1);
            output = combustiontoolbox.gui.AppController.evaluateConsoleCommands(obj.app, commands);
        end

        function setup = currentProblemSetup(obj)
            if ~isempty(obj.app)
                input = obj.readInput();
                currentMixture = [];

                if isprop(obj.app, 'mixture')
                    currentMixture = obj.app.mixture;
                end

                problem = obj.solver.buildProblem(input, currentMixture);
                setup = problem.setup;
                return
            end

            if isempty(obj.currentSolution) || ~isfield(obj.currentSolution, 'problem') ...
                    || ~isfield(obj.currentSolution.problem, 'setup')
                error('AppController:MissingProblemSetup', ...
                    'Solve a problem before exporting a MATLAB script.');
            end

            setup = obj.currentSolution.problem.setup;
        end

        function filepath = selectProblemScriptFile(~, setup)
            problemType = matlab.lang.makeValidName(lower(char(setup.problemType)));
            defaultName = ['ct_', problemType, '_setup.m'];
            filter = {'*.m', 'MATLAB script (*.m)'};
            [filename, folder] = uiputfile(filter, 'Export problem setup', defaultName);

            if isequal(filename, 0)
                filepath = [];
                return
            end

            filepath = fullfile(folder, filename);
        end

        function filepath = ensureMatlabScriptExtension(~, filepath)
            [folder, name, extension] = fileparts(filepath);

            if isempty(extension)
                filepath = fullfile(folder, [name, '.m']);
            end
        end

        function writeTextFile(~, filepath, text)
            [fileId, message] = fopen(filepath, 'w');

            if fileId < 0
                error('AppController:ScriptExportFailed', ...
                    'Could not write MATLAB script "%s": %s', filepath, message);
            end

            cleanup = onCleanup(@() fclose(fileId));
            fwrite(fileId, text, 'char');
            clear cleanup
        end

        function setScriptExportStatus(obj, filepath)
            if isempty(obj.app) || isempty(filepath)
                return
            end

            obj.setStatus('done');
            obj.consolePanel.setOutput(sprintf('Problem script exported to:\n%s', filepath));
        end

        function output = formatConsoleError(~, exception)
            stackName = 'unknown';
            stackLine = 0;

            if ~isempty(exception.stack)
                stackName = exception.stack(1).name;
                stackLine = exception.stack(1).line;
            end

            output = sprintf('Error in function %s() at line %d.\nIdentifier: %s\nMessage: %s', ...
                stackName, stackLine, exception.identifier, exception.message);
        end
    end

    methods (Static, Access = private)
        function output = evaluateConsoleCommands(app, commands) %#ok<INUSD>
            output = '';

            for i = 1:numel(commands)
                output = evalc(char(commands{i}));
            end
        end
    end

    methods (Access = private)
        function setStatus(obj, status)
            if isempty(obj.app)
                return
            end

            if isempty(obj.statusPanel)
                return
            end

            switch status
                case 'idle'
                    obj.statusPanel.setIdle();
                case 'working'
                    obj.statusPanel.setWorking();
                case 'done'
                    obj.statusPanel.setDone();
                case 'error'
                    obj.statusPanel.setError();
            end
        end

        function setProblemErrorStatus(obj)
            if ~isempty(obj.statusPanel)
                obj.statusPanel.setProblemError();
            end
        end

        function clearProblemErrorStatus(obj)
            if ~isempty(obj.statusPanel)
                obj.statusPanel.clearProblemError();
            end
        end

        function setFinishedStatus(obj)
            if isempty(obj.app)
                return
            end

            problemError = obj.problemErrorValue();

            if problemError > obj.app.maxRelativeError
                obj.setStatus('error');
                obj.setProblemErrorStatus();
                obj.consolePanel.setOutput(sprintf( ...
                    'Warning! The maximum relative error is %.2f%%. Results may be compromised.\nDecreasing the tolerance and increasing the number of iterations may solve the problem.', ...
                    problemError * 100));
                return
            end

            obj.setStatus('done');
            obj.clearProblemErrorStatus();
            obj.consolePanel.setOutput('Done! check tab "Results"');
        end

        function value = problemErrorValue(obj)
            value = 0;

            if ~isempty(obj.resultsPanel)
                value = obj.resultsPanel.problemErrorValue();
            end
        end

        function setErrorStatus(obj, exception)
            if isempty(obj.app)
                rethrow(exception);
            end

            obj.setStatus('error');
            [message, title] = obj.errorStatusMessage(exception);
            obj.consolePanel.setOutput(message);

            if obj.shouldShowAlert()
                uialert(obj.app.UIFigure, {message}, title, 'Icon', 'warning');
            end
        end

        function [message, title] = errorStatusMessage(~, exception)
            if strcmp(exception.identifier, 'AppProblemBuilder:MissingReactants')
                message = exception.message;
                title = 'Setup required';
                return
            end

            stackName = 'unknown';
            stackLine = 0;

            if ~isempty(exception.stack)
                stackName = exception.stack(1).name;
                stackLine = exception.stack(1).line;
            end

            message = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                stackName, stackLine, exception.message);
            title = 'Warning';
        end

        function value = shouldShowAlert(obj)
            value = isprop(obj.app, 'UIFigure') && ~isempty(obj.app.UIFigure) ...
                && isvalid(obj.app.UIFigure) ...
                && strcmpi(char(obj.app.UIFigure.Visible), 'on');
        end
    end

end

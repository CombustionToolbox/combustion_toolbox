classdef AppConsolePanel < handle
    % Manages the GUI command window controls
    %
    % AppConsolePanel owns command-window presentation details so the main
    % App Designer class can delegate output, history, and geometry updates
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %     commandHistory (cell): MATLAB command history entries
    %     commandHistoryIndex (double): Current command history index
    %
    % Example:
    %     * panel = combustiontoolbox.gui.AppConsolePanel(app);

    properties (Access = private)
        app
        commandHistory
        commandHistoryIndex
    end

    methods
        function obj = AppConsolePanel(app)
            % AppConsolePanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppConsolePanel): Initialized command window panel
            if nargin < 1
                app = [];
            end

            obj.app = app;
            obj.commandHistory = {};
            obj.commandHistoryIndex = 0;
        end

        function initialize(obj, welcomeText)
            % Initialize output and command history
            %
            % Args:
            %     welcomeText (char): Initial command-window output
            obj.setOutput(welcomeText);
            obj.refreshHistory();
        end

        function restoreOutput(obj, previousValue)
            % Restore output when the user tries to edit it
            %
            % Args:
            %     previousValue (cell): Previous command-window output
            obj.app.Console_text.Value = previousValue;
        end

        function maximize(obj)
            % Maximize the command-window output area
            obj.app.Console_text.Position = [78, 30, 554, 462];
        end

        function minimize(obj)
            % Minimize the command-window output area
            obj.app.Console_text.Position = [78, 30, 554, 57];
        end

        function clearOutput(obj)
            % Clear command-window output
            obj.app.Console_text.Value = '';
        end

        function setOutput(obj, value)
            % Set command-window output
            %
            % Args:
            %     value (char | string | cell): Output text
            obj.app.Console_text.Value = value;
        end

        function selectHistory(obj, key)
            % Select a command-history entry from a keyboard event
            %
            % Args:
            %     key (char): Key identifier from UIFigureKeyPress
            switch lower(key)
                case 'uparrow'
                    obj.commandHistoryIndex = obj.commandHistoryIndex - 1;
                case 'downarrow'
                    obj.commandHistoryIndex = obj.commandHistoryIndex + 1;
                case 'rightarrow'
                    obj.refreshHistory();
                case 'leftarrow'
                    obj.commandHistoryIndex = 1;
            end

            obj.applyHistorySelection();
        end

        function refreshHistory(obj)
            % Refresh MATLAB command history
            obj.commandHistory = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
            obj.commandHistoryIndex = length(obj.commandHistory);
        end
    end

    methods (Access = private)
        function applyHistorySelection(obj)
            % Apply the selected history entry to the input field
            if obj.commandHistoryIndex >= 1 && obj.commandHistoryIndex <= numel(obj.commandHistory)
                obj.app.Console.Value = char(obj.commandHistory(obj.commandHistoryIndex));
                return
            end

            obj.app.Console.Value = 'Value out of current command history.';
            obj.refreshHistory();
        end
    end
end

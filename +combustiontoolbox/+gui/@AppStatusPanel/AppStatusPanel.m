classdef AppStatusPanel < handle
    % Applies GUI status colors and problem-error highlighting
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %
    % Examples:
    %     * panel = combustiontoolbox.gui.AppStatusPanel(app);
    %     * panel.setWorking();

    properties (Constant, Access = private)
        idleColor        = [0.8000, 0.8000, 0.8000]
        workingColor     = [0.9961, 0.9804, 0.8314]
        doneColor        = [0.5608, 0.7255, 0.6588]
        errorColor       = [0.6400, 0.0800, 0.1800]
        warningColor     = [0.9961, 0.9804, 0.8314]
        defaultTextColor = [0, 0, 0]
    end

    properties (Access = private)
        app
    end

    methods
        function obj = AppStatusPanel(app)
            % AppStatusPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppStatusPanel): Initialized status panel
            if nargin < 1
                app = [];
            end

            obj.app = app;
        end

        function setIdle(obj)
            % Set idle status
            obj.setLampColor(obj.idleColor);
        end

        function setWorking(obj)
            % Set working status
            obj.setLampColor(obj.workingColor);
        end

        function setDone(obj)
            % Set done status
            obj.setLampColor(obj.doneColor);
        end

        function setError(obj)
            % Set error status
            obj.setLampColor(obj.errorColor);
        end

        function setWarning(obj)
            % Set warning status
            obj.setLampColor(obj.warningColor);
        end

        function setProblemError(obj)
            % Highlight problem result error status
            obj.setColor('text_error_problem', 'FontColor', obj.errorColor);
            obj.setColor('ResultsTab', 'ForegroundColor', obj.errorColor);
        end

        function clearProblemError(obj)
            % Restore problem result error status
            obj.setColor('text_error_problem', 'FontColor', obj.defaultTextColor);
            obj.setColor('ResultsTab', 'ForegroundColor', obj.defaultTextColor);
        end
    end

    methods (Access = private)
        function setLampColor(obj, color)
            % Set lamp color when the component exists
            if obj.hasComponent('Lamp') && isprop(obj.app.Lamp, 'Color')
                obj.app.Lamp.Color = color;
            end
        end

        function setColor(obj, componentName, propertyName, color)
            % Set component color when the component and property exist
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), propertyName)
                obj.app.(componentName).(propertyName) = color;
            end
        end

        function value = hasComponent(obj, componentName)
            % Check whether an App Designer component is available
            value = isobject(obj.app) && isprop(obj.app, componentName) ...
                && ~isempty(obj.app.(componentName));
        end
    end
end

classdef AppProblemInputsPanel < handle
    % Coordinates coupled problem input controls
    %
    % AppProblemInputsPanel owns small interactions between PR and PP inputs,
    % including reactant/product pressure mirroring, shock Mach/velocity
    % coupling, and mutually exclusive fields for rocket, detonation, and
    % shock variants
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %     mirroredPressureProblems (cell): Problem types sharing pressure input
    %
    % Examples:
    %     * panel = combustiontoolbox.gui.AppProblemInputsPanel(app);
    %     * panel.onReactantMachChanged();

    properties
        app
    end

    properties (Access = private, Constant)
        mirroredPressureProblems = {'TP', 'HP'}
    end

    methods
        function obj = AppProblemInputsPanel(app)
            % AppProblemInputsPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppProblemInputsPanel): Initialized problem inputs panel
            if nargin < 1
                app = [];
            end

            obj.app = app;
        end

        function onProductPressureChanged(obj, event)
            % Mirror product pressure into reactant pressure when linked
            if obj.shouldMirrorReactantProductPressure()
                obj.setValue('PR2', obj.eventValue(event, obj.componentValue('PP2', '')));
            end
        end

        function onReactantPressureChanged(obj, event)
            % Mirror reactant pressure and refresh shock velocity
            if obj.shouldMirrorReactantProductPressure()
                obj.setValue('PP2', obj.eventValue(event, obj.componentValue('PR2', '')));
            end

            obj.updateMachOrVelocity('Mach');
        end

        function onReactantMachChanged(obj)
            % Refresh shock velocity from Mach number
            obj.updateMachOrVelocity('Mach');
        end

        function onReactantVelocityChanged(obj)
            % Refresh shock Mach number from velocity
            obj.updateMachOrVelocity('velocity');
        end

        function onRocketReactantAreaChanging(obj)
            % Keep rocket area ratio inputs mutually exclusive
            obj.clearIfProblem('ROCKET', 'PP3');
        end

        function onRocketProductAreaChanging(obj)
            % Keep rocket area ratio inputs mutually exclusive
            obj.clearIfProblem('ROCKET', 'PR3');
        end

        function onDetonationReactantAngleChanging(obj)
            % Keep detonation angle inputs mutually exclusive
            obj.clearIfAnyProblem({'DET_OBLIQUE', 'DET_POLAR'}, 'PP4');
        end

        function onDetonationProductAngleChanging(obj)
            % Keep detonation angle inputs mutually exclusive
            obj.clearIfAnyProblem({'DET_OBLIQUE', 'DET_POLAR'}, 'PR4');
        end

        function onShockReactantAngleChanging(obj)
            % Keep shock angle inputs mutually exclusive
            obj.clearIfAnyProblem({'SHOCK_OBLIQUE', 'SHOCK_POLAR'}, 'PP5');
        end

        function onShockProductAngleChanging(obj)
            % Keep shock angle inputs mutually exclusive
            obj.clearIfAnyProblem({'SHOCK_OBLIQUE', 'SHOCK_POLAR'}, 'PR5');
        end

        function updateMachOrVelocity(obj, inputName)
            % Update shock velocity from Mach number or Mach number from velocity
            %
            % Args:
            %     inputName (char): Source input name, `Mach` or `velocity`
            if ~obj.problemContains('SHOCK')
                return
            end

            if isempty(obj.componentData('UITable_R', {}))
                obj.setValue('PR3', '-');
                return
            end

            mix = obj.app.mixture;

            switch lower(inputName)
                case 'mach'
                    [mach, FLAG_ARRAY] = combustiontoolbox.gui.AppInput.parseValue(obj.componentValue('PR4', []));
                    velocity = mach * mix.sound;
                    obj.setValue('PR3', obj.formatVectorOrScalar(velocity, FLAG_ARRAY));
                otherwise
                    [velocity, FLAG_ARRAY] = combustiontoolbox.gui.AppInput.parseValue(obj.componentValue('PR3', []));
                    mach = velocity / mix.sound;
                    obj.setValue('PR4', obj.formatVectorOrScalar(mach, FLAG_ARRAY));
            end
        end
    end

    methods (Access = private)
        function value = shouldMirrorReactantProductPressure(obj)
            problemType = char(obj.componentValue('ProblemType', ''));
            value = any(strcmpi(problemType, obj.mirroredPressureProblems));
        end

        function clearIfProblem(obj, problemType, componentName)
            if strcmpi(char(obj.componentValue('ProblemType', '')), problemType)
                obj.setValue(componentName, '');
            end
        end

        function clearIfAnyProblem(obj, problemTypes, componentName)
            for i = 1:numel(problemTypes)
                if obj.problemContains(problemTypes{i})
                    obj.setValue(componentName, '');
                    return
                end
            end
        end

        function value = problemContains(obj, pattern)
            value = contains(char(obj.componentValue('ProblemType', '')), pattern, 'IgnoreCase', true);
        end

        function value = componentValue(obj, componentName, defaultValue)
            value = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.app.(componentName).Value;
            end
        end

        function data = componentData(obj, componentName, defaultValue)
            data = defaultValue;

            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Data')
                data = obj.app.(componentName).Data;
            end
        end

        function setValue(obj, componentName, value)
            if obj.hasComponent(componentName) && isprop(obj.app.(componentName), 'Value')
                value = obj.assignmentValue(obj.app.(componentName), value);
                obj.app.(componentName).Value = value;
            end
        end

        function value = eventValue(obj, event, defaultValue) %#ok<INUSD>
            value = defaultValue;

            if isempty(event)
                return
            end

            if isstruct(event) && isfield(event, 'Value')
                value = event.Value;
                return
            end

            if isobject(event) && isprop(event, 'Value')
                value = event.Value;
            end
        end

        function value = formatVectorOrScalar(obj, value, FLAG_ARRAY) %#ok<INUSD>
            if FLAG_ARRAY
                value = sprintf('%.2f:%.2f:%.2f', value(1), value(2) - value(1), value(end));
            else
                value = sprintf('%.2f', value);
            end
        end

        function value = hasComponent(obj, componentName)
            value = isobject(obj.app) && isprop(obj.app, componentName) && ~isempty(obj.app.(componentName));
        end

        function value = assignmentValue(~, component, value)
            if ~isnumeric(component.Value)
                return
            end

            if isempty(value)
                value = 0;
            elseif ischar(value) || isstring(value)
                numericValue = str2double(value);

                if isnan(numericValue)
                    value = 0;
                else
                    value = numericValue;
                end
            elseif ~isnumeric(value)
                value = 0;
            end

            if isnumeric(value) && ~isempty(value)
                value = value(1);
            end
        end
    end
end

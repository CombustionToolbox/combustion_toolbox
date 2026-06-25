classdef AppSettingsPanel < handle
    % Applies quick settings controls to GUI services
    %
    % AppSettingsPanel owns the link between App Designer settings controls
    % and CT solver or plotting services
    %
    % Attributes:
    %     app (combustion_toolbox): Main App Designer object
    %     session (AppSession): Long-lived GUI services
    %
    % Example:
    %     * panel = combustiontoolbox.gui.AppSettingsPanel(app, session);

    properties (Access = private)
        app
        session
    end

    methods
        function obj = AppSettingsPanel(app, session)
            % AppSettingsPanel constructor
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %     session (AppSession): Long-lived GUI services
            %
            % Returns:
            %     obj (AppSettingsPanel): Initialized settings panel
            if nargin < 1
                app = [];
            end

            if nargin < 2
                session = [];
            end

            obj.app = app;
            obj.session = session;
        end

        function applyTraceTolerance(obj)
            % Apply equilibrium trace tolerance
            obj.session.equilibriumSolver.tolMoles = obj.value('TraceoptionEditField');
        end

        function applyRootFindingTolerance(obj)
            % Apply equilibrium root-finding tolerance
            obj.session.equilibriumSolver.tol0 = obj.value('RootFindingMethodEditField');
        end

        function applyShockDetonationTolerance(obj)
            % Apply shock and detonation tolerance
            value = obj.value('ShocksandDetonationsEditField');
            obj.session.shockSolver.tol0 = value;
            obj.session.detonationSolver.tol0 = value;
        end

        function applyDisplaySpeciesTolerance(obj)
            % Apply minimum displayed species tolerance
            obj.session.plotConfig.mintolDisplay = obj.value('DisplaySpeciesEditField');
        end

        function applyEquilibriumMaxIterations(obj)
            % Apply equilibrium maximum iterations
            obj.session.equilibriumSolver.itMax = obj.value('MaxiterationsRFMEditField');
        end

        function applyShockDetonationMaxIterations(obj)
            % Apply shock and detonation maximum iterations
            value = obj.value('MaxiterationsSDEditField');
            obj.session.shockSolver.itMax = value;
            obj.session.detonationSolver.itMax = value;
        end

        function applyRootTemperatureLeft(obj)
            % Apply left root temperature bound
            obj.session.equilibriumSolver.root_T0_l = obj.value('RFMT0_LEditField');
        end

        function applyRootTemperatureRight(obj)
            % Apply right root temperature bound
            obj.session.equilibriumSolver.root_T0_r = obj.value('RFMT0_REditField');
        end

        function applyRootTemperatureInitial(obj)
            % Apply initial root temperature
            obj.session.equilibriumSolver.root_T0 = obj.value('RFMT0EditField');
        end

        function refreshFromPreference(obj, tag)
            % Refresh quick settings controls affected by a preferences tag
            if isempty(obj.session)
                return
            end

            tag = obj.normalizeTag(tag);

            if numel(tag) < 2
                return
            end

            switch strjoin(tag(1:2), '/')
                case 'equilibriumSolver/tolMoles'
                    obj.setControlValue('TraceoptionEditField', ...
                        obj.session.equilibriumSolver.tolMoles);
                case 'equilibriumSolver/tol0'
                    obj.setControlValue('RootFindingMethodEditField', ...
                        obj.session.equilibriumSolver.tol0);
                case 'equilibriumSolver/itMax'
                    obj.setControlValue('MaxiterationsRFMEditField', ...
                        obj.session.equilibriumSolver.itMax);
                case 'equilibriumSolver/root_T0_l'
                    obj.setControlValue('RFMT0_LEditField', ...
                        obj.session.equilibriumSolver.root_T0_l);
                case 'equilibriumSolver/root_T0_r'
                    obj.setControlValue('RFMT0_REditField', ...
                        obj.session.equilibriumSolver.root_T0_r);
                case 'equilibriumSolver/root_T0'
                    obj.setControlValue('RFMT0EditField', ...
                        obj.session.equilibriumSolver.root_T0);
                case 'plotConfig/mintolDisplay'
                    obj.setControlValue('DisplaySpeciesEditField', ...
                        obj.session.plotConfig.mintolDisplay);
                case {'shockSolver/tol0', 'detonationSolver/tol0'}
                    obj.refreshSharedShockDetonationControl('ShocksandDetonationsEditField', 'tol0');
                case {'shockSolver/itMax', 'detonationSolver/itMax'}
                    obj.refreshSharedShockDetonationControl('MaxiterationsSDEditField', 'itMax');
            end
        end

        function refreshAllFromSession(obj)
            % Refresh all quick settings controls from session services
            if isempty(obj.session)
                return
            end

            obj.setControlValue('TraceoptionEditField', ...
                obj.session.equilibriumSolver.tolMoles);
            obj.setControlValue('RootFindingMethodEditField', ...
                obj.session.equilibriumSolver.tol0);
            obj.setControlValue('DisplaySpeciesEditField', ...
                obj.session.plotConfig.mintolDisplay);
            obj.setControlValue('MaxiterationsRFMEditField', ...
                obj.session.equilibriumSolver.itMax);
            obj.setControlValue('RFMT0_LEditField', ...
                obj.session.equilibriumSolver.root_T0_l);
            obj.setControlValue('RFMT0_REditField', ...
                obj.session.equilibriumSolver.root_T0_r);
            obj.setControlValue('RFMT0EditField', ...
                obj.session.equilibriumSolver.root_T0);
            obj.refreshSharedShockDetonationControl('ShocksandDetonationsEditField', 'tol0');
            obj.refreshSharedShockDetonationControl('MaxiterationsSDEditField', 'itMax');
        end
    end

    methods (Access = private)
        function value = value(obj, componentName)
            % Return the value of a settings component
            value = obj.app.(componentName).Value;
        end

        function setControlValue(obj, componentName, value)
            % Set a quick settings value when the control exists
            if isempty(obj.app) || ~isprop(obj.app, componentName) ...
                    || isempty(obj.app.(componentName)) ...
                    || ~isprop(obj.app.(componentName), 'Value')
                return
            end

            obj.app.(componentName).Value = value;
        end

        function refreshSharedShockDetonationControl(obj, componentName, propertyName)
            % Refresh a shared shock and detonation setting when values agree
            shockValue = obj.session.shockSolver.(propertyName);
            detonationValue = obj.session.detonationSolver.(propertyName);

            if isequal(shockValue, detonationValue)
                obj.setControlValue(componentName, shockValue);
            end
        end

        function tag = normalizeTag(~, tag)
            % Normalize a preference tag into a row cell array
            if isstring(tag)
                tag = cellstr(tag(:))';
            elseif ischar(tag)
                tag = strsplit(tag, filesep);
            elseif iscell(tag)
                tag = tag(:)';
            else
                tag = {};
            end
        end
    end
end

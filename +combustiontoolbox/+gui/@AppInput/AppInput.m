classdef AppInput < handle
    % Collects, normalizes, and parses Combustion Toolbox GUI inputs.
    %
    % Attributes:
    %     problemType (char): Selected problem type identifier
    %     reactants: Current reactants dropdown value
    %     products: Current products dropdown value
    %     productSpecies (cell): Species considered as products
    %     displaySpecies (cell): Species selected for composition plots
    %     displaySpeciesTolerance (float): Minimum species fraction shown in plots
    %     reportType (char): Report generation mode selected in the GUI
    %     reactantsTable (cell): Data from the reactants UITable
    %     temperatureReactants (float): Reactants temperature input [K]
    %     pressureReactants (float): Reactants pressure input [bar]
    %     temperatureProducts (float): Products temperature input [K]
    %     pressureProducts (float): Products pressure input [bar]
    %     equivalenceRatio (float): Equivalence ratio input [-]
    %     flags (struct): Boolean GUI options
    %     additionalProperties (struct): Problem-specific PR and PP properties
    %
    % Examples:
    %     * input = combustiontoolbox.gui.AppInput.fromApp(app);
    %     * value = combustiontoolbox.gui.AppInput.parseValue('300:100:700');

    properties
        problemType = 'TP'
        reactants
        products
        productSpecies = {}
        displaySpecies = {}
        displaySpeciesTolerance = 1e-6
        reportType = 'Auto'
        reactantsTable = {}
        temperatureReactants = []
        pressureReactants = []
        temperatureProducts = []
        pressureProducts = []
        equivalenceRatio = []
        flags = struct()
        additionalProperties = struct()
    end

    properties (Access = private)
        rawValues = struct()
    end

    methods
        function obj = AppInput(varargin)
            % AppInput constructor
            %
            % Optional Args:
            %     data (struct): Struct with fields matching AppInput properties
            %
            % Returns:
            %     obj (AppInput): Initialized input snapshot
            if nargin == 0
                return
            end

            data = varargin{1};

            if ~isstruct(data)
                error('AppInput:InvalidInput', 'AppInput constructor expects a struct.');
            end

            fields = fieldnames(data);

            for i = 1:numel(fields)

                if ~isprop(obj, fields{i})
                    continue
                end

                obj.(fields{i}) = data.(fields{i});
            end

        end

        function expressions = exportExpressions(obj)
            % Return raw expressions that should be preserved in scripts
            %
            % Returns:
            %     expressions (struct): Numeric values and raw expressions
            names = fieldnames(obj.rawValues);
            expressions = struct('value', {}, 'text', {});
            count = 0;

            for i = 1:numel(names)
                rawValue = obj.rawValues.(names{i});
                parsedValue = combustiontoolbox.gui.AppInput.parseValue(rawValue);
                rawExpression = combustiontoolbox.gui.AppInput.rawExpression(rawValue);

                if isempty(rawExpression) || numel(parsedValue) < 2
                    continue
                end

                count = count + 1;
                expressions(count).value = parsedValue;
                expressions(count).text = rawExpression;
            end
        end
    end

    methods (Static)
        function obj = fromApp(app)
            % Build a typed input snapshot from the current App Designer app
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppInput): Typed snapshot of the current GUI inputs
            data.problemType = combustiontoolbox.gui.AppInput.getComponentValue(app, 'ProblemType', 'TP');
            data.reactants = combustiontoolbox.gui.AppInput.getDropdownText(app, 'Reactants', []);
            data.products = combustiontoolbox.gui.AppInput.getDropdownText(app, 'Products', []);
            data.productSpecies = combustiontoolbox.gui.AppInput.getComponentItems(app, 'listbox_Products');
            data.displaySpecies = combustiontoolbox.gui.AppInput.getComponentItems(app, 'listbox_LS_display');
            data.displaySpeciesTolerance = combustiontoolbox.gui.AppInput.getComponentValue(app, 'DisplaySpeciesEditField', 1e-6);
            data.reportType = combustiontoolbox.gui.AppInput.getComponentValue(app, 'Report_type', 'Auto');
            data.reactantsTable = combustiontoolbox.gui.AppInput.getComponentData(app, 'UITable_R');
            data.rawValues = combustiontoolbox.gui.AppInput.rawComponentValues(app);

            data.temperatureReactants = combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PR1');
            data.pressureReactants = combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PR2');
            data.temperatureProducts = combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP1');
            data.pressureProducts = combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP2');
            data.equivalenceRatio = combustiontoolbox.gui.AppInput.parseComponentValue(app, 'edit_phi');
            data.additionalProperties = struct( ...
                'PR3', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PR3'), ...
                'PR4', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PR4'), ...
                'PR5', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PR5'), ...
                'PP3', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP3'), ...
                'PP4', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP4'), ...
                'PP5', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP5'), ...
                'PP6', combustiontoolbox.gui.AppInput.parseComponentValue(app, 'PP6'));

            data.flags = struct( ...
                'frozenChemistry', combustiontoolbox.gui.AppInput.getComponentValue(app, 'FrozenchemistryCheckBox', false), ...
                'ionizedSpecies', combustiontoolbox.gui.AppInput.getComponentValue(app, 'IonizedspeciesCheckBox', false), ...
                'idealAir', combustiontoolbox.gui.AppInput.getComponentValue(app, 'IdealAirCheckBox', false), ...
                'infiniteAreaChamber', combustiontoolbox.gui.AppInput.getComponentValue(app, 'FLAG_IAC', true), ...
                'printResults', combustiontoolbox.gui.AppInput.getComponentValue(app, 'PrintresultsCheckBox', false));

            obj = combustiontoolbox.gui.AppInput(data);
        end

        function [value, FLAG_ARRAY] = parseValue(rawValue, varargin)
            % Parse a GUI text value into numeric data without eval/str2num
            %
            % Supported forms:
            %   300
            %   [300 500 700]
            %   [300, 500, 700]
            %   300:100:700
            %
            % Optional selectors:
            %   'first', 'last'
            %   'number', N, 'direction', 'last'
            %
            % Args:
            %     rawValue (char | string | numeric | cell): GUI value to parse
            %
            % Optional Args:
            %     direction (char): 'first' or 'last'
            %     number (float): Number of values to return
            %
            % Returns:
            %     value (float): Parsed numeric value
            %     FLAG_ARRAY (logical): True when rawValue represents an array
            FLAG_ARRAY = false;

            if nargin == 0 || isempty(rawValue)
                value = [];
                return
            end

            if iscell(rawValue)
                if isempty(rawValue)
                    value = [];
                    return
                end

                rawValue = rawValue{1};
            end

            if isnumeric(rawValue)
                value = rawValue;
                FLAG_ARRAY = numel(value) > 1;
            else
                textValue = strtrim(char(rawValue));

                if isempty(textValue) || strcmp(textValue, '-')
                    value = [];
                    return
                end

                [value, FLAG_ARRAY] = combustiontoolbox.gui.AppInput.parseTextValue(textValue);
            end

            if isempty(varargin) || isempty(value)
                return
            end

            [n, direction] = combustiontoolbox.gui.AppInput.parseSelectionOptions(varargin{:});
            value = combustiontoolbox.gui.AppInput.selectValues(value, n, direction);
        end
    end

    methods (Static, Access = private)
        function value = getComponentValue(app, componentName, defaultValue)
            value = defaultValue;
            [component, exists] = combustiontoolbox.gui.AppInput.getComponent(app, componentName);

            if exists && combustiontoolbox.gui.AppInput.hasComponentProperty(component, 'Value')
                value = component.Value;
            end
        end

        function items = getComponentItems(app, componentName)
            items = {};
            [component, exists] = combustiontoolbox.gui.AppInput.getComponent(app, componentName);

            if exists && combustiontoolbox.gui.AppInput.hasComponentProperty(component, 'Items')
                items = component.Items;
            end
        end

        function value = getDropdownText(app, componentName, defaultValue)
            value = combustiontoolbox.gui.AppInput.getComponentValue(app, componentName, defaultValue);
            [component, exists] = combustiontoolbox.gui.AppInput.getComponent(app, componentName);

            if ~exists || ~combustiontoolbox.gui.AppInput.hasComponentProperty(component, 'Items')
                return
            end

            items = component.Items;

            if isempty(items)
                return
            end

            if combustiontoolbox.gui.AppInput.hasComponentProperty(component, 'ItemsData') ...
                    && ~isempty(component.ItemsData)
                index = combustiontoolbox.gui.AppInput.dropdownItemsDataIndex( ...
                    value, component.ItemsData);

                if ~isempty(index) && index <= numel(items)
                    if iscell(items)
                        value = items{index};
                    else
                        value = items(index);
                    end

                    if isstring(value)
                        value = char(value);
                    end

                    return
                end
            end

            if ischar(value) || isstring(value)
                value = char(value);
            end
        end

        function index = dropdownItemsDataIndex(value, itemsData)
            index = [];

            if isnumeric(itemsData) || islogical(itemsData)
                if ~(isnumeric(value) || islogical(value))
                    return
                end

                index = find(itemsData == value, 1);
                return
            end

            if iscell(itemsData)
                index = find(cellfun(@(item) isequal(item, value), itemsData), 1);
                return
            end

            if isstring(itemsData)
                if ~(ischar(value) || isstring(value))
                    return
                end

                index = find(strcmp(cellstr(itemsData), char(value)), 1);
            end
        end

        function data = getComponentData(app, componentName)
            data = {};
            [component, exists] = combustiontoolbox.gui.AppInput.getComponent(app, componentName);

            if exists && combustiontoolbox.gui.AppInput.hasComponentProperty(component, 'Data')
                data = component.Data;
            end
        end

        function [component, exists] = getComponent(app, componentName)
            component = [];
            exists = isobject(app) && isprop(app, componentName);

            if exists
                component = app.(componentName);
            end
        end

        function value = hasComponentProperty(component, propertyName)
            value = isobject(component) && isprop(component, propertyName);
        end

        function value = parseComponentValue(app, componentName)
            rawValue = combustiontoolbox.gui.AppInput.getComponentValue(app, componentName, []);
            value = combustiontoolbox.gui.AppInput.parseValue(rawValue);
        end

        function values = rawComponentValues(app)
            componentNames = {'PR1', 'PR2', 'PR3', 'PR4', 'PR5', ...
                'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'PP6', 'edit_phi'};
            values = struct();

            for i = 1:numel(componentNames)
                name = componentNames{i};
                values.(name) = combustiontoolbox.gui.AppInput.getComponentValue(app, name, []);
            end
        end

        function [value, FLAG_ARRAY] = parseTextValue(textValue)
            FLAG_ARRAY = false;

            if startsWith(textValue, '[') && endsWith(textValue, ']')
                FLAG_ARRAY = true;
                textValue = extractBetween(string(textValue), 2, strlength(textValue) - 1);
                textValue = char(textValue);

                if count(string(textValue), ":") == 2
                    values = combustiontoolbox.gui.AppInput.extractNumbers(textValue);

                    if numel(values) == 3
                        value = values(1):values(2):values(3);
                        return
                    end
                end

                value = combustiontoolbox.gui.AppInput.extractNumbers(textValue);
                return
            end

            if count(string(textValue), ":") == 2
                values = combustiontoolbox.gui.AppInput.extractNumbers(textValue);

                if numel(values) == 3
                    FLAG_ARRAY = true;
                    value = values(1):values(2):values(3);
                    return
                end
            end

            value = combustiontoolbox.gui.AppInput.extractNumbers(textValue);
        end

        function numbers = extractNumbers(textValue)
            numberPattern = '[-+]?(?:(?:\d+\.?\d*)|(?:\.\d+))(?:[eEdD][-+]?\d+)?';
            tokens = regexp(textValue, numberPattern, 'match');

            if isempty(tokens)
                numbers = [];
                return
            end

            tokens = strrep(tokens, 'D', 'E');
            tokens = strrep(tokens, 'd', 'e');
            numbers = str2double(tokens);
        end

        function expression = rawExpression(rawValue)
            expression = '';

            if isstring(rawValue) && isscalar(rawValue)
                rawValue = char(rawValue);
            end

            if ~ischar(rawValue)
                return
            end

            expression = strtrim(rawValue);
        end

        function [n, direction] = parseSelectionOptions(varargin)
            n = 1;
            direction = 'first';
            i = 1;

            while i <= numel(varargin)
                option = varargin{i};

                if ~(ischar(option) || isstring(option))
                    i = i + 1;
                    continue
                end

                switch lower(char(option))
                    case {'first', 'last'}
                        direction = lower(char(option));
                        i = i + 1;
                    case {'number', 'n', 'nindex'}
                        n = varargin{i + 1};
                        i = i + 2;
                    case {'direction', 'get', 'position', 'pos'}
                        direction = lower(char(varargin{i + 1}));
                        i = i + 2;
                    otherwise
                        i = i + 1;
                end
            end
        end

        function value = selectValues(value, n, direction)
            value = value(:).';

            if isempty(value)
                return
            end

            n = min(n, numel(value));

            switch lower(direction)
                case 'last'
                    value = value(end - n + 1:end);
                otherwise
                    value = value(1:n);
            end
        end
    end
end

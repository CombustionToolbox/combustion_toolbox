classdef AppState < handle
    % Stores mutable GUI state independently of App Designer handles.
    %
    % Attributes:
    %     selectedProblem (char): Selected problem type identifier
    %     selectedReactants: Current reactants dropdown value
    %     selectedProducts: Current products dropdown value
    %     productSpecies (cell): Species considered as products
    %     displaySpecies (cell): Species selected for composition plots
    %     flags (struct): Boolean GUI options
    %     hiddenTabs (struct): Tabs hidden by the GUI
    %     lastResultsId: Identifier of the last selected result set
    %
    % Examples:
    %     * state = combustiontoolbox.gui.AppState(input);
    %     * state.updateFromInput(input);

    properties
        selectedProblem = 'TP'
        selectedReactants = []
        selectedProducts = []
        productSpecies = {}
        displaySpecies = {}
        flags = struct()
        hiddenTabs = struct()
        lastResultsId = []
    end

    methods
        function obj = AppState(varargin)
            % AppState constructor
            %
            % Optional Args:
            %     input (AppInput | struct): Initial values for the state
            %
            % Returns:
            %     obj (AppState): Initialized state object
            if nargin == 0
                return
            end

            input = varargin{1};

            if isa(input, 'combustiontoolbox.gui.AppInput')
                updateFromInput(obj, input);
            elseif isstruct(input)
                updateFromStruct(obj, input);
            else
                error('AppState:InvalidInput', 'AppState expects an AppInput or struct.');
            end
        end

        function updateFromInput(obj, input)
            % Update state from a typed AppInput snapshot
            %
            % Args:
            %     input (AppInput): Current GUI input snapshot
            obj.selectedProblem = input.problemType;
            obj.selectedReactants = input.reactants;
            obj.selectedProducts = input.products;
            obj.productSpecies = input.productSpecies;
            obj.displaySpecies = input.displaySpecies;
            obj.flags = input.flags;
        end

        function updateFromStruct(obj, input)
            % Update matching state properties from a struct
            %
            % Args:
            %     input (struct): Struct with fields matching AppState properties
            fields = fieldnames(input);

            for i = 1:numel(fields)
                if isprop(obj, fields{i})
                    obj.(fields{i}) = input.(fields{i});
                end
            end
        end
    end
end

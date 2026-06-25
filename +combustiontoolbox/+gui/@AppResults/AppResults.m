classdef AppResults < handle
    % GUI-side storage for computed studies and selected result views
    %
    % AppResults keeps the current result set separate from the App Designer
    % tree, tables, axes, and text fields. This makes it possible to test result
    % selection, output labeling, and serialization without creating UI components
    %
    % Attributes:
    %     items (struct): Current GUI-friendly results
    %     temporary (struct): Last raw temporary results array
    %     layout (struct): Result layout metadata used by the current results
    %     selectedIndex (float): Selected result index
    %
    % Examples:
    %     * results = combustiontoolbox.gui.AppResults();
    %     * results.set(resultArray);
    %     * selected = results.selected();

    properties
        items = struct([])
        temporary = struct([])
        layout = struct()
        selectedIndex = 1
    end

    methods
        function set(obj, results, varargin)
            % Store a new result set
            %
            % Args:
            %     results (struct): GUI-friendly result array
            %
            % Optional Args:
            %     temporary (struct): Raw temporary result array
            obj.items = results;

            if nargin > 2
                obj.temporary = varargin{1};
            else
                obj.temporary = results;
            end

            obj.selectedIndex = 1;
        end

        function setFromSolverOutputs(obj, problemType, solverOutputs, layout, varargin)
            % Store ordered solver outputs using result layout metadata
            %
            % Args:
            %     problemType (char): Problem type identifier
            %     solverOutputs (cell): Ordered solver outputs
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %
            % Optional Args:
            %     variant (char): Active output variant
            results = obj.buildFromSolverOutputs(problemType, solverOutputs, layout, varargin{:});
            obj.layout = layout;
            obj.set(results, results);
        end

        function results = buildFromSolverOutputs(obj, problemType, solverOutputs, layout, varargin)
            % Build GUI result structs from ordered solver outputs
            %
            % Args:
            %     problemType (char): Problem type identifier
            %     solverOutputs (cell): Ordered solver outputs
            %     layout (struct): Result layout metadata from AppProblemCatalog
            %
            % Optional Args:
            %     variant (char): Active output variant
            %
            % Returns:
            %     results (struct): Per-case GUI result data
            variant = '';

            if nargin > 4
                variant = varargin{1};
            end

            if ~iscell(solverOutputs)
                solverOutputs = {solverOutputs};
            end

            outputStates = obj.selectOutputStates(layout, variant);

            if numel(solverOutputs) < numel(outputStates)
                error('AppResults:MissingSolverOutput', ...
                    'Expected %d solver outputs for %s, received %d.', ...
                    numel(outputStates), problemType, numel(solverOutputs));
            end

            numCases = obj.numberOfCases(solverOutputs, numel(outputStates));
            results = repmat(struct(), 1, numCases);

            for i = numCases:-1:1
                results(i).ProblemType = problemType;
                results(i).outputStates = outputStates;
                results(i).states = struct();

                for j = 1:numel(outputStates)
                    outputValue = obj.valueForCase(solverOutputs{j}, i);
                    results(i).(outputStates(j).field) = outputValue;
                    results(i).states.(outputStates(j).id) = outputValue;
                end
            end
        end

        function clear(obj)
            % Clear all stored results and reset selection
            obj.items = struct([]);
            obj.temporary = struct([]);
            obj.layout = struct();
            obj.selectedIndex = 1;
        end

        function value = hasResults(obj)
            % Check if the object contains results
            %
            % Returns:
            %     value (logical): True if at least one result is stored
            value = ~isempty(obj.items);
        end

        function result = selected(obj)
            % Return the selected result
            %
            % Returns:
            %     result (struct): Selected result, or [] if no result exists
            if ~hasResults(obj)
                result = [];
                return
            end

            index = min(obj.selectedIndex, numel(obj.items));
            result = obj.items(index);
        end
    end

    methods (Access = private)
        function outputStates = selectOutputStates(~, layout, variant)
            outputStates = layout.outputStates;

            if isempty(variant) || ~isfield(layout, 'variantOutputStates')
                return
            end

            variant = char(variant);
            variantNames = fieldnames(layout.variantOutputStates);

            for i = 1:numel(variantNames)
                variantName = variantNames{i};

                if strcmpi(variant, variantName) || endsWith(upper(variant), ['_', upper(variantName)])
                    outputStates = layout.variantOutputStates.(variantName);
                    return
                end
            end
        end

        function numCases = numberOfCases(~, solverOutputs, numOutputStates)
            counts = zeros(1, numOutputStates);

            for i = 1:numOutputStates
                counts(i) = numel(solverOutputs{i});
            end

            numCases = max(counts);
        end

        function value = valueForCase(~, value, index)
            if isempty(value) || isscalar(value)
                return
            end

            value = value(min(index, numel(value)));
        end
    end
end

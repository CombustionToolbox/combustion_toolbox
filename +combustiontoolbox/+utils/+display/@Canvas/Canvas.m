classdef (Abstract) Canvas < handle
    % The :mat:class:`Canvas` class is an abstract base class used to manage 
    % figure creation and plotting layout configuration in the :mat:pkg:`+display` module.
    %
    % The :mat:class:`Canvas` class standardizes figure and axis formatting, 
    % provides tiled layout management, and offers utilities for axis selection and cycling.
    %
    % A :mat:class:`Canvas` subclass can be initialized as follows: ::
    %
    %      config = combustiontoolbox.utils.display.PlotConfig();
    %      fig = combustiontoolbox.utils.display.PlotFigure(config);
    %
    % See also: :mat:class:`PlotConfig`, :mat:class:`PlotFigure`, :mat:class:`PlotComposition`

    properties
        ax matlab.graphics.axis.Axes                             % Axis handle
        tiled matlab.graphics.layout.TiledChartLayout            % Tiled layout for multiple axes
        fig matlab.ui.Figure                                     % Figure handle
        config (1, 1) combustiontoolbox.utils.display.PlotConfig % Plot configuration object
    end

    properties (Access = private)
        currentIndex = 1 % Current index for sequential filling of axes
    end

    properties (Dependent)
        size    % Size of the axis array
        numAxis % Number of axes in the figure
    end

    methods 
    
        function value = get.size(obj)
            % Get the size axis array
            value = size(obj.ax);
        end

        function value = get.numAxis(obj)
            % Get the number of axes in the figure
            value = numel(obj.ax);
        end

    end

    methods (Access = public)
        
        function obj = Canvas(config)
            % Constructor
            if nargin > 0
                obj.config = config;
                return
            end
            
            % Default configuration
            obj.config = combustiontoolbox.utils.display.PlotConfig();
        end

        function obj = add(obj)
            % Add a new axis to the figure
            % This method is intended to be overridden in subclasses
            % to add specific axes or plots.
    
            % Definitions
            config = obj.config;

            % Check if fig is already created
            if isempty(obj.fig)
                obj.fig = figure;
                set(obj.fig, 'units', 'normalized', ...
                    'innerposition', config.innerpositionLayout, ...
                    'outerposition', config.outerpositionLayout)
            end

            % Check if tiled layout is already created
            if isempty(obj.tiled)
                obj.tiled = tiledlayout(obj.fig, 'flow', ...
                            'Padding', obj.config.padding, ...
                            'TileSpacing', obj.config.tilespacing);
            end
            
            % Create a new axis in the tiled layout
            ax = nexttile(obj.tiled);

            % Set the axes properties
            setFormat(obj, ax);

            % Append new axis
            obj.ax(end + 1) = ax;

            % Add a listener so when the user closes the figure, we reset
            addlistener(obj.fig, 'ObjectBeingDestroyed', @(~,~) obj.onFigureClose());
        end

        function [ax, config, fig] = new(obj)
            % Create and configure a new figure

            % Definitions
            config = obj.config;

            % Create figure
            fig = figure;
            set(fig, 'units', 'normalized', ...
                'innerposition', config.innerposition, ...
                'outerposition', config.outerposition)
                
            % Create axes
            ax = axes(fig);

            % Set axes properties
            setFormat(obj, ax);

            % Store the axes and figure in the object
            obj.ax = ax;
            obj.fig = fig;

            % Add a listener so when the user closes the figure, we reset
            addlistener(fig, 'ObjectBeingDestroyed', @(~,~) obj.onFigureClose());
        end

        function [axArray, config, fig] = newArray(obj, varargin)
            % Create a figure with an nrows x ncols tiled layout
            % and return an array of PlotFigure subclass instances

            % Definitions
            config = obj.config;

            % Initialization
            nrows = [];
            ncols = [];

            % Parse input arguments
            if nargin > 1
                nrows = varargin{1};
                ncols = varargin{2};
            end

            % Create parent figure and tiled layout
            fig = figure;
            set(fig, 'units', 'normalized', ...
                'innerposition', config.innerpositionLayout, ...
                'outerposition', config.outerpositionLayout)

            % Check if nrows and ncols are provided
            if isempty(nrows) || isempty(ncols)
                nrows = 1; ncols = 1;
                tiled = tiledlayout(fig, 'flow', ...
                                'Padding', config.padding, ...
                                'TileSpacing', config.tilespacing);
            else
                tiled = tiledlayout(fig, nrows, ncols, ...
                                'Padding', config.padding, ...
                                'TileSpacing', config.tilespacing);
            end

            % Initialize array of objects
            axArray = gobjects(ncols, nrows); % Preallocate array of axes

            % Loop through rows and columns to create axes
            for i = 1:ncols
                
                for j = 1:nrows
                    % Create new axes in the tiled layout
                    ax = nexttile(tiled);

                    % Set the axes properties
                    setFormat(obj, ax);

                    % Save the axis in the array
                    axArray(i, j) = ax;
                end

            end

            % Store the axes and figure in the object
            obj.ax = axArray;  
            obj.tiled = tiled;
            obj.fig = fig;

            % Add a listener so when the user closes the figure, we reset
            addlistener(fig, 'ObjectBeingDestroyed', @(~,~) obj.onFigureClose());
        end

        function next(obj)
            % Cycle to the next index
            obj.currentIndex = mod(obj.currentIndex, obj.numAxis) + 1;
        end

        function ax = getAxis(obj, varargin)
            % Return axis to use for plotting
            %
            % Args:
            %     obj (PlotFigure): Instance of the PlotFigure superclass
            %
            % Optional Args:
            %     index (vector): Index of the axis to use [row, col]
            %     FLAG_SEQUENTIAL (bool): If true, use sequential filling of axes
            %
            % Returns:
            %     ax (Axes): Axis to use for plotting
            %
            % Examples:
            %     * ax = obj.getAxis();
            %     * ax = obj.getAxis('index', [1, 2]);
            %     * ax = obj.getAxis('FLAG_SEQUENTIAL', true);
            %     * ax = obj.getAxis('index', [1, 2], 'FLAG_SEQUENTIAL', true); % Start sequence at (1, 2)

            % Check if obj.ax exists and is valid
            if isempty(obj.ax)
                ax = new(obj);
                return
            end

            % If varargin is empty, return current axis
            if nargin < 1
                ax = obj.ax(1, 1);
                return
            end

            % Parse input
            p = inputParser;
            addParameter(p, 'index', [], @(x) isvector(x) && length(x) == 2);
            addParameter(p, 'FLAG_SEQUENTIAL', false, @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            % Set properties
            index = p.Results.index;
            FLAG_SEQUENTIAL = p.Results.FLAG_SEQUENTIAL;

            % Get axis based on index
            if ~FLAG_SEQUENTIAL
                if isempty(index)
                    index1 = ceil(obj.currentIndex / obj.size(2));
                    index2 = mod(obj.currentIndex - 1, obj.size(2)) + 1;
                end

                ax = obj.ax(index1, index2);
                return
            end

            % Use sequential filling if multiple axes exist
            if ~isempty(index)
                obj.currentIndex = sub2ind(size(obj.ax), index(1), index(2));
            end

            % Get axis index
            k = obj.currentIndex;

            % Get axis
            ax = obj.ax(k);

            % Increment the current index for sequential filling
            next(obj)
        end

        function ax = setFormat(obj, ax)
            % Set the format of the axes based on the configuration
            % 
            % Args:
            %     obj (PlotFigure): Instance of the PlotFigure superclass
            %     ax (Axes): Axes to format
            %
            % Returns:
            %     ax (Axes): Formatted axes
            %
            % Example:
            %     ax = obj.setFormat(gca);

            % Definitions
            config = obj.config;

            % Set axes properties
            set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize - 2);
            set(ax, 'BoxStyle', 'full', 'TickLabelInterpreter', 'latex');
            set(ax, 'Layer', 'Top');
            set(ax, 'XScale', config.xscale, 'YScale', config.yscale);
            set(ax, 'XDir', config.xdir, 'YDir', config.ydir);
            xlim(ax, config.axis_x);
            ylim(ax, config.axis_y);
            box(ax, config.box);
            grid(ax, config.grid);
            hold(ax, config.hold);
            xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'Interpreter', 'latex');
            ylabel(ax, config.labely, 'FontSize', config.fontsize, 'Interpreter', 'latex');
        end

        function setTitle(obj, titleText)
            % Set figure title with LaTeX interpreter
            if nargin < 2 || isempty(titleText)
                titleText = obj.config.title;
            end

            if isempty(titleText)
                return
            end

            title(obj.ax, titleText, 'Interpreter', 'latex', 'FontSize', obj.config.fontsize + 4);
        end

        function setLegend(obj, legendLabels)
            % Set legend on the axis
            legend(obj.ax, legendLabels, 'FontSize', obj.config.fontsize - 2, ...
                'Location', obj.config.legend_location, 'Interpreter', 'latex');
        end
        
    end

    methods (Abstract)
        plot(obj, varargin)
    end

    methods (Access = private)

        function onFigureClose(obj)
            % Reset the axes and figure when the figure is closed
            obj.ax = matlab.graphics.axis.Axes.empty(1,0);
            obj.tiled = matlab.graphics.layout.TiledChartLayout.empty(1,0);
            obj.fig = matlab.ui.Figure.empty(1,0);
            obj.currentIndex = 1;
        end
    
    end

end
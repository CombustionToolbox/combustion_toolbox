classdef PlotFigure < combustiontoolbox.utils.display.Canvas
    % The :mat:func:`PlotFigure` class extends :mat:func:`Canvas` to provide a unified
    % plotting interface for thermodynamic and transport properties (e.g., `cp`, `h`, `Xi`)
    % against an arbitrary independent variable.
    %
    % A :mat:func:`PlotFigure` object can be initialized as follows: ::
    %
    %      config = combustiontoolbox.utils.display.PlotConfig();
    %      fig = combustiontoolbox.utils.display.PlotFigure(config);
    %      fig.plot('equivalenceRatio', mixArray, 'T', mixArray);
    %
    % This creates a formatted plot with consistent styling using :mat:func:`PlotConfig`.
    %
    % See also: :mat:func:`Canvas`, :mat:func:`PlotComposition`, :mat:func:`PlotConfig`

    properties (Access = private)
        colorIndex = 1; % Index for color cycling
    end

    methods

        function obj = PlotFigure(config)
            % Constructor
            if nargin < 1
                config = combustiontoolbox.utils.display.PlotConfig();
            end

            obj@combustiontoolbox.utils.display.Canvas(config);
        end

        function plot(obj, x_field, x_var, y_field, y_var, varargin)
            % Plot method for scalar or multiple y_fields
            %
            % Args:
            %     x_field (char): Field name for x-axis
            %     x_var (cell): Cell array of structures/objects with x_field
            %     y_fields (char or cell): Field name(s) for y-axis
            %     y_var (cell): Cell array with same length as x_var
            %     varargin: Optional key-value arguments
            %
            % Supported options:
            %     'config', 'legend', 'legend_location', 'linestyle',
            %     'linewidth', 'fontsize', 'title', 'xlabel', 'ylabel',
            %     'label_type', 'xscale', 'yscale', 'xdir', 'ydir',
            %     'color', 'basis', 'validation', 'nfrec'
    
            % Import packages
            import combustiontoolbox.utils.display.*

            % Default settings
            FLAG_BASIS = false;
            FLAG_COLOR_NEW = false;
            config = obj.config;
            ax = obj.ax;
        
            % Get x and y values
            try
                x = cell2vector(x_var, x_field);
                y = cell2vector(y_var, y_field);
            catch
                x = x_var.x_field;
                y = y_var.y_field;
            end
        
            % Check lenghts
            if length(x) < 2
                return
            end
        
            % Check aditional inputs
            for i = 1:2:nargin - 5
        
                switch lower(varargin{i})
                    case 'config'
                        config = varargin{i + 1};
                        config.labelx = interpreterLabel(x_field, config.label_type, false);
                        config.labely = interpreterLabel(y_field, config.label_type, false);
                    case {'leg', 'legend'}
                        config.legend_name = varargin{i + 1};
                    case {'legend_location'}
                        config.legend_location = varargin{i + 1};
                    case {'ax', 'axes'}
                        ax = varargin{i + 1};
                    case 'linestyle'
                        config.linestyle = varargin{i + 1};
                    case 'linewidth'
                        config.linewidth = varargin{i + 1};
                    case 'fontsize'
                        config.fontsize = varargin{i + 1};
                    case 'title'
                        config.title = varargin{i + 1};
                    case {'labelx', 'xlabel', 'label_x', 'x_label'}
                        config.labelx = varargin{i + 1};
                    case {'labely', 'ylabel', 'label_y', 'y_label'}
                        config.labelx = varargin{i + 1};
                    case {'label_type'}
                        config.label_type = varargin{i + 1};
                    case {'xscale'}
                        config.xscale = varargin{i + 1};
                    case {'yscale'}
                        config.yscale = varargin{i + 1};
                    case {'xdir'}
                        config.xdir = varargin{i + 1};
                    case {'ydir'}
                        config.ydir = varargin{i + 1};
                    case 'color'
                        config.colorline = varargin{i + 1};
        
                        if ~isfloat(config.colorline)
                            FLAG_COLOR_NEW = true;
                        end
                    case 'basis'
                        basis = varargin{i + 1};
        
                        if ~isempty(basis)
                            FLAG_BASIS = true;
                        end
                        
                end
        
            end
            
            if isempty(config.labelx)
                config.labelx = interpreterLabel(x_field, config.label_type, false);
            end
        
            if isempty(config.labely)
                config.labely = interpreterLabel(y_field, config.label_type, false);
            end
        
            % Create figure (if necessary)
            if isempty(ax)
                obj.config = config;
                ax = obj.new();
            else
                ax = getAxis(obj);
            end
            
            % change units if required
            switch y_field
                case {'cp', 'cv', 'hf', 'ef', 'h', 'e', 'g', 's'}
                    y = y * 1e-3; % [kJ ...]
            end
        
            % Check if property has to be divided by the basis (kg or mol)
            if FLAG_BASIS
                y_basis = cell2vector(y_var, basis);
                y = y ./ y_basis;
        
                config.labelx = interpreterLabel(x_field, config.label_type, true, basis);
        
                if ~strcmpi(config.labely, 'Multiple variables')
                    config.labely = interpreterLabel(y_field, config.label_type, true, basis);
                end
        
            end
        
            % Plot
            if FLAG_COLOR_NEW
                dline = plot(ax, x, y, config.linestyle, 'LineWidth', config.linewidth);
            else
                dline = plot(ax, x, y, config.linestyle, 'LineWidth', config.linewidth, 'Color', config.colorline);
            end
            
            % Set labels
            xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'interpreter', 'latex');
            ylabel(ax, config.labely, 'FontSize', config.fontsize, 'interpreter', 'latex');

            % Title and legend
            if ~isempty(config.title)
                setTitle(obj, config.title);
            end

            if ~isempty(config.legend_name)
                setLegend(obj, config.legend_name);
            end

        end

        function plotArray(obj, x_field, x_var, y_field, y_var, varargin)
            % Plot method for multiple y_fields in an array format
            %
            % Args:
            %     x_field (char): Field name for x-axis
            %     x_var (cell): Cell array of structures/objects with x_field
            %     y_fields (char or cell): Field name(s) for y-axis
            %     y_var (cell): Cell array with same length as x_var
            %     varargin: Optional key-value arguments
            
            % Import packages
            import combustiontoolbox.utils.display.*

            % Check inputs 
            if ~iscell(x_field), x_field = {x_field}; end
            if ~iscell(y_field), y_field = {y_field}; end

            % Default settings
            nfrec = 1;
            FLAG_PLOT_VALIDATION = false;
            FLAG_BASIS = false;
            FLAG_SAME = false;
            FLAG_COLOR = false;
            tiled = obj.tiled;
            config = obj.config;
            
            % Definitions
            cmap = config.cmap;
            selectColor = obj.colorIndex;
            selectSymbol = 1;
            
            % Check aditional inputs
            for i = 1:2:nargin - 5

                switch lower(varargin{i})
                    case {'validation', 'results'}
                        results2 = varargin{i + 1};
                        FLAG_PLOT_VALIDATION = true;
                    case 'config'
                        config = varargin{i + 1};
                    case {'leg', 'legend'}
                        config.legend_name = varargin{i + 1};
                    case {'legend_location'}
                        config.legend_location = varargin{i + 1};
                    case {'ax', 'axes'}
                        tiled = varargin{i + 1};
                    case 'linestyle'
                        config.linestyle = varargin{i + 1};
                    case 'linewidth'
                        config.linewidth = varargin{i + 1};
                    case 'fontsize'
                        config.fontsize = varargin{i + 1};
                    case 'title'
                        config.title = varargin{i + 1};
                    case {'labelx', 'xlabel', 'label_x', 'x_label'}
                        config.labelx = varargin{i + 1};
                    case {'labely', 'ylabel', 'label_y', 'y_label'}
                        config.labelx = varargin{i + 1};
                    case {'label_type'}
                        config.label_type = varargin{i + 1};
                    case {'xscale'}
                        config.xscale = varargin{i + 1};
                    case {'yscale'}
                        config.yscale = varargin{i + 1};
                    case {'xdir'}
                        config.xdir = varargin{i + 1};
                    case {'ydir'}
                        config.ydir = varargin{i + 1};
                    case 'color'
                        config.colorline = varargin{i + 1};
                        cmap = config.colorline;
                        FLAG_COLOR = true;
                    case 'basis'
                        basis = varargin{i + 1};
                        FLAG_BASIS = true;
                    case 'nfrec'
                        nfrec = varargin{i + 1};
                end

            end

            % Create main figure
            if isempty(tiled)
                obj.newArray();
            else
                FLAG_SAME = true;
                
                if ~FLAG_COLOR
                    nextColor(obj);
                    selectColor = obj.colorIndex;
                end

            end

            % Set common title
            if isa(x_var, 'combustiontoolbox.core.Mixture')
                if isempty(config.title)
                    config.title = getTitle(x_var(1));
                else
                    config.title = sprintf('%s - %s', config.title, getTitle(x_var(1)));
                end
            end

            setTitle(tiled, config)
            config.title = [];

            % Definitions
            N_properties = length(y_field);
            if ~FLAG_BASIS
                basis = cell(1, N_properties);
            end

            % Plot properties
            for i = 1:N_properties
                
                if FLAG_SAME
                    obj.next();
                elseif i > 1
                    obj.add(); obj.next();
                end

                ax = getAxis(obj);
                obj.plot(x_field{i}, x_var, y_field{i}, y_var, 'config', config, 'ax', ax, 'basis', basis{i}, 'color', cmap(selectColor, :));

                if FLAG_PLOT_VALIDATION
                    plot(ax, results2.(x_field{i})(1:nfrec:end), results2.(y_field{i})(1:nfrec:end), config.symbolStyles{selectSymbol}, 'LineWidth', config.linewidth, 'color', 'k', 'MarkerFaceColor', cmap(selectColor, :));
                end

            end

        end

        function nextColor(obj)
            % Cycle to the next color in the color palette
            obj.colorIndex = mod(obj.colorIndex, size(obj.config.cmap, 1)) + 1;
        end

    end

end
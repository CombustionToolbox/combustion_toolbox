function [main_ax, main_figure] = plotProperties(x_field, x_var, y_field, y_var, varargin)
    % Plot figure with customizable settings
    %
    % Args:
    %     x_field (cell): Cell array containing field names for the x-axis data
    %     x_var (cell): Cell array containing the x-axis data
    %     y_field (cell): Cell array containing field names for the y-axis data
    %     y_var (cell): Cell array containing the y-axis data
    %
    % Optional Name-Value Pair Args:
    %     * config (struct): Struct with the configuration for plots
    %     * leg or legend (cell): Cell array of strings containing the legend names
    %     * legend_location (char): Location of the legend
    %     * ax or axes (object): Handle of the axes to plot on
    %     * linestyle (char): Line style
    %     * linewidth (float): Line width
    %     * fontsize (float): Font size
    %     * title (char): Title of the figure
    %     * labelx, xlabel, label_x, or x_label (char): x-axis label
    %     * labely, ylabel, label_y, or 'y_label (char): y-axis label
    %     * label_type (char): Label type
    %     * xscale (char): Set x-axis scale (linear or log)
    %     * yscale (char): Set y-axis scale (linear or log)
    %     * xdir (char): Set x-axis direction (normal or reverse)
    %     * ydir (char): Set y-axis direction (normal or reverse)
    %     * color (float): Line color [R, G, B]
    %
    % Returns:
    %     mainfigure (obj): Figure object
    
    % Import packages
    import combustiontoolbox.utils.display.*

    % Default settings
    nfrec = 1;
    FLAG_PLOT_VALIDATION = false;
    FLAG_BASIS = false;
    FLAG_SAME = false;
    main_ax = [];
    config = PlotConfig();

    % Definitions
    colorPalette = config.colorlines(3, :);
    lineStyles = {'-', '--', ':', '-.'};
    symbolStyles = {'d', 'o', 's', '<'};
    selectColor = 1;
    selectLine = 1;
    selectSymbol = 1;

    % Check inputs 
    if ~iscell(x_field), x_field = {x_field}; end
    if ~iscell(y_field), y_field = {y_field}; end

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
                main_ax = varargin{i + 1};
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
            case 'basis'
                basis = varargin{i + 1};
                FLAG_BASIS = true;
            case 'nfrec'
                nfrec = varargin{i + 1};
        end

    end
    
    % Create main figure
    if isempty(main_ax)
        main_figure = figure;
        set(main_figure, 'units', 'normalized', 'innerposition', config.innerposition, 'outerposition', config.outerposition);
        main_ax = tiledlayout(main_figure, 'flow');
    else
        main_figure = gcf;
        FLAG_SAME = true;
    end

    % Set common title
    setTitle(main_ax, config)
    config.title = [];
    
    % Definitions
    N_properties = length(y_field);
    if ~FLAG_BASIS
        basis = cell(1, N_properties);
    end

    % Plot properties
    for i = 1:N_properties
        
        if FLAG_SAME
            nexttile(main_ax, i)
            ax = gca;
            setFigure(ax, config);
        else
            nexttile;
            ax = gca;
            setFigure(ax, config);
        end

        plotFigure(x_field{i}, x_var, y_field{i}, y_var, 'config', config, 'ax', ax, 'basis', basis{i}, 'color', 'color', colorPalette(selectColor, :));

        if FLAG_PLOT_VALIDATION
            plot(ax, results2.(x_field{i})(1:nfrec:end), results2.(y_field{i})(1:nfrec:end), symbolStyles{selectSymbol}, 'LineWidth', config.linewidth, 'color', colorPalette(selectColor, :), 'MarkerFaceColor', 'white');
        end

    end
    

end
function ax = plot_figure(x_field, x_var, y_field, y_var, varargin)
    % Plot figure

    % Default settings
    FLAG_COLOR_NEW = false;
    Misc = Miscellaneous();
    config = Misc.config;
    config.labelx = interpreter_label(x_field, config.label_type); % Set x label
    config.labely = interpreter_label(y_field, config.label_type); % Set y label
    ax = []; % Set empty axes
    % Get x and y values
    x = cell2vector(x_var, x_field);
    y = cell2vector(y_var, y_field);
    % Check aditional inputs
    for i = 1:2:nargin - 5

        switch lower(varargin{i})
            case 'config'
                config = varargin{i + 1};
                config.labelx = interpreter_label(x_field, config.label_type);
                config.labely = interpreter_label(y_field, config.label_type);
            case {'leg', 'legend'}
                config.legend_name = varargin{i + 1};
            case {'legend_location'}
                config.legend_location = varargin{i + 1};
            case {'ax', 'axes', 'figure'}
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

        end

    end

    % Create figure (if necessary)
    if isempty(ax)
        ax = set_figure(config);
    end

    % Plot
    if FLAG_COLOR_NEW
        plot(ax, x, y, config.linestyle, 'LineWidth', config.linewidth);
    else
        plot(ax, x, y, config.linestyle, 'LineWidth', config.linewidth, 'Color', config.colorline);
    end
    
    % Set labels
    xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'interpreter', 'latex');
    ylabel(ax, config.labely, 'FontSize', config.fontsize, 'interpreter', 'latex');

    % Set title
    if ~isempty(config.title)
        title(ax, config.title, 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
    end

    % Legend
    if ~isempty(config.legend_name)
        legend(ax, config.legend_name, 'FontSize', config.fontsize - 2, 'Location', config.legend_location, 'interpreter', 'latex');
    end

end

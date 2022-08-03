function ax = plot_figure(x_field, x_var, y_field, y_var, varargin)
    % Plot figure

    % Default settings
    Misc = Miscellaneous();
    config = Misc.config;
    config.labelx = interpreter_label(x_field, config.label_type); % Set x label
    config.labely = interpreter_label(y_field, config.label_type); % Set y label
    ax = []; % Set empty axes 
    % Get x and y values
    x = cell2vector(x_var, x_field);
    y = cell2vector(y_var, y_field);
    % Check aditional inputs
    for i=1:2:nargin-5
        switch lower(varargin{i})
            case 'config'
                config = varargin{i+1};
                if isempty(config.labelx)
                    config.labelx = interpreter_label(x_field, config.label_type);
                end
                if isempty(config.labely)
                    config.labely = interpreter_label(y_field, config.label_type);
                end
            case {'leg', 'legend'}
                config.legend_name = varargin{i+1};
            case {'ax', 'axes', 'figure'}
                ax = varargin{i+1};
            case 'linewidth'
                config.linewidth = varargin{i+1};
            case 'fontsize'
                config.fontsize = varargin{i+1};
            case 'title'
                config.title = varargin{i+1};
            case {'labelx', 'xlabel', 'label_x', 'x_label'}
                config.labelx = varargin{i+1};
            case {'labely', 'ylabel', 'label_y', 'y_label'}
                config.labelx = varargin{i+1};
            case {'label_type'}
                config.label_type = varargin{i+1};

        end
    end
    % Create figure (if necessary)
    if isempty(ax)
        ax = set_figure(config);
    end
    % Plot
    plot(ax, x, y, 'LineWidth', config.linewidth);
    % Legend
    if ~isempty(config.legend_name)
        legend(ax, config.legend_name, 'FontSize', config.fontsize-2, 'Location', 'northeastoutside', 'interpreter', 'latex');
    end
end
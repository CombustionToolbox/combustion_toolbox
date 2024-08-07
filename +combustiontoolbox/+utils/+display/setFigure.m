function [ax, config, fig] = setFigure(varargin)
    % Initialize figure with a standard composition
    %
    % Optional Args:
    %     * ax (axis): Figure axis
    %     * config (struct): Struct with default plot parameters
    %
    % Returns:
    %    Tuple containing
    %
    %     * ax (axis): Axis of the standard figure
    %     * config (struct): Struct with default plot parameters
    %     * fig (figure): Standard figure
    %
    % Example:
    %     [ax, config, fig] = setFigure();
    %     [ax, config, fig] = setFigure(ax);
    %     [ax, config, fig] = setFigure(config);
    %     [ax, config, fig] = setFigure(ax, config);
    %
    % See also: plot_settings
    
    % Import packages
    import combustiontoolbox.utils.display.PlotConfig

    % Default values
    config = PlotConfig();
    ax = [];
    
    % Unpack input
    if nargin > 0

        for i = 1:nargin

            if isa(varargin{i}, 'matlab.graphics.axis.Axes')
                ax = varargin{i};
            elseif isa(varargin{i}, 'combustiontoolbox.utils.display.PlotConfig')
                config = varargin{i};
            end

        end

        if ~isempty(ax)

            if isempty(ax.XLabel)
                config.labelx = ax.XLabel.String;
            end

            if isempty(ax.YLabel)
                config.labely = ax.YLabel.String;
            end

        end

    end

    % Set axes
    if isempty(ax)
        fig = figure;
        set(fig, 'units', 'normalized', 'innerposition', config.innerposition, ...
            'outerposition', config.outerposition)
        ax = axes(fig);
    end

    set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize - 2)
    set(ax, 'BoxStyle', 'full', 'TickLabelInterpreter', 'latex')
    set(ax, 'Layer', 'Top');
    set(ax, 'xscale', config.xscale, 'yscale', config.yscale)
    set(ax, 'XDir', config.xdir, 'YDir', config.ydir)
    xlim(ax, config.axis_x);
    ylim(ax, config.axis_y);
    box(ax, config.box);
    grid(ax, config.grid);
    hold(ax, config.hold);
    xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'interpreter', 'latex');
    ylabel(ax, config.labely, 'FontSize', config.fontsize, 'interpreter', 'latex');
end
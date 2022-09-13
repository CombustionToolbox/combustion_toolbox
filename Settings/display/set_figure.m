function [ax, config, fig] = set_figure(varargin)
    % Initialize figure with a standard composition
    %
    % Optional args:
    %   * ax (axis):       Figure axis
    %   * config (struct): Struct with default plot parameters
    %
    % Returns:
    %   ax (axis):         Axis of the standard figure
    %   config (struct):   Struct with default plot parameters
    %   fig (figure):      Standard figure

    % Default values
    Misc = Miscellaneous();
    config = Misc.config;
    ax = [];
    % Unpack input
    if nargin > 0
        for i=1:nargin
            if isobject(varargin{i})
                ax = varargin{i};
            else
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
    if isempty(ax)
        fig = figure;
        set(fig,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
            'outerposition',[0.1 0.1 0.9 0.8])
        ax = axes(fig);
    end
    
    set(ax,'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    xlim(ax, config.axis_x);
    ylim(ax, config.axis_y);
    box(ax, config.box);
    grid(ax, config.grid);
    hold(ax, config.hold);
    xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'interpreter', 'latex');
    ylabel(ax, config.labely, 'FontSize', config.fontsize, 'interpreter', 'latex');
%     title({strcat('$',config.tit,'$')},'Interpreter','latex','FontSize',config.fontsize+4)
end
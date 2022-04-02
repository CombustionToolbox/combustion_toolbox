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
    % Unpack input
    if nargin > 0
        ax = varargin{1};
        if nargin > 1
            config = varargin{2};
        end
        if ~isempty(ax.XLabel)
            config.labelx = ax.XLabel;
        end
        if ~isempty(ax.YLabel)
            config.labely = ax.YLabel;
        end
    else
        fig = figure;
        set(fig,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
            'outerposition',[0.1 0.1 0.9 0.8])
        ax = axes(fig);
    end
    
    set(ax,'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
    hold(ax, 'on'); axis(ax, 'tight');
    xlabel(ax, config.labelx, 'FontSize', config.fontsize, 'interpreter', 'latex');
    ylabel(ax, config.labely, 'FontSize', config.fontsize, 'interpreter', 'latex');
%     title({strcat('$',config.tit,'$')},'Interpreter','latex','FontSize',config.fontsize+4)
end
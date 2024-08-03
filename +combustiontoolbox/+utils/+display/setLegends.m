function setLegends(ax, legend_name, varargin)
    % Set legend to the given axes
    %
    % Args:
    %     ax (object): Handle to the axes
    %     legend_name (cell): Cell array of char containing the legend names
    %
    % Optional Name-Value Pairs Args:
    %     * config (struct): Struct containing the configuration parameters for the plots
    %     * obj (object): Handle to the plotted objects (e.g. lines, patches, etc.)
    
    % Import packages
    import combustiontoolbox.utils.display.PlotConfig

    % Default values
    config = PlotConfig();
    FLAG_OBJECTS = false;
    
    % Unpack inputs
    for i = 1:2:nargin-2
        switch lower(varargin{i})
            case {'config'}
                config = varargin{i + 1};
            case {'objects', 'obj'}
                obj = varargin{i + 1};
        end
    end

    if FLAG_OBJECTS
        legend(ax, obj, legend_name, 'FontSize', config.fontsize - 4, 'Location', 'best', 'interpreter', 'latex');
        return
    end
    
    legend(ax, legend_name, 'FontSize', config.fontsize - 4, 'Location', 'best', 'interpreter', 'latex');
end

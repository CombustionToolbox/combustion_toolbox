function set_legends(ax, legend_name, varargin)
    % Set legend to the given axes
    if nargin < 3
        config.fontsize = 18;
    else
        config = varargin{1};
    end

    legend_name = strrep(legend_name, '_', '\_');
    legend(ax, legend_name, 'FontSize', config.fontsize - 4, 'Location', 'best', 'interpreter', 'latex');
end

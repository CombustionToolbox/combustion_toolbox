function set_legends(ax, legend_name, varargin)
    % Set legend to the given axes
    if nargin < 3
        config.fontsize  = 18;
    else
        config = varargin{1};
    end
    legend(ax, legend_name,'FontSize', config.fontsize-2, 'Location', 'best', 'interpreter', 'latex');
end
function set_title(ax, varargin)
    % Set legend to the given axes
    if nargin < 2
        config.fontsize  = 18;
    else
        config = varargin{1};
    end
    title(ax, config.tit, 'Interpreter', 'latex', 'FontSize', config.fontsize+4);
end
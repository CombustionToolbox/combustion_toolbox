function set_legends(ax, legend_name, varargin)
    % Set legend to the given axes

    % Default values
    config.fontsize = 18;
    FLAG_SPECIES = true;
    % Unpack inputs
    for i = 1:2:nargin-2
        switch lower(varargin{i})
            case {'config'}
                config = varargin{i + 1};
            case {'flag', 'flag_species'}
                FLAG_SPECIES = varargin{i + 1};
        end
    end

    if FLAG_SPECIES
        legend_name = strrep(legend_name, '_', '\_');
    end

    legend(ax, legend_name, 'FontSize', config.fontsize - 4, 'Location', 'best', 'interpreter', 'latex');
end

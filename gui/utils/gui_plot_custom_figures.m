function gui_plot_custom_figures(app)
    % Function that plot custom figures based on the parameters of the
    % GUI's custom figures tab
    
    % Clear axes
    cla(app.UIAxes);
    % Get mixtures
    mixtures_nodes = app.Tree_mixtures.CheckedNodes;
    % Get variables x
    x_field = get(app.Tree_variable_x.CheckedNodes, 'Text');
    % Get variables y
    y_field = get(app.Tree_variable_y.CheckedNodes, 'Text');
    % Check non empty values
    if isempty(x_field) || isempty(y_field), return; end
    % Check that variables y is a cell
    if ~iscell(y_field), y_field = {y_field}; end
    % Plot settings
    Misc = Miscellaneous();
    config = Misc.config;
    config.labelx = interpreter_label(x_field, config.label_type); % Set x label
    config.labely = interpreter_label(y_field, config.label_type); % Set y label
    NUM_NODES = length(mixtures_nodes);
    NUM_PROP = length(y_field);
    LINE_STYLES = {'-', '--', ':', '-.'};
    if NUM_PROP > 1
        COLOR_PALETTE = brewermap(NUM_PROP, config.colorpalette);
    else
        COLOR_PALETTE = config.colorline;
    end
    legend_name = cell(1, NUM_NODES * NUM_PROP);
    legend_name_1 = get(mixtures_nodes, 'Text');
    legend_name_2 = y_field;
    % Check that legend_name_1 is a cell
    if ~iscell(legend_name_1), legend_name_1 = {legend_name_1}; end
    % Initialize axes
    if app.DefaultsettingsCheckBox.Value
        ax = set_figure(app.UIAxes, config);
    end
    % Plot
    k = NUM_NODES * NUM_PROP;
    for i = NUM_NODES:-1:1
        mixtures = app.temp_results.(mixtures_nodes(i).Text);
        x_values = cell2vector(mixtures, x_field);
        for j = NUM_PROP:-1:1
            plot_figure(x_field, x_values, y_field{j}, mixtures, 'config', config, 'axes', app.UIAxes, 'linestyle', LINE_STYLES{i}, 'color', COLOR_PALETTE(j, :));
            legend_name{k} = [legend_name_1{i}, ' - ', legend_name_2{j}];
            k = k - 1;
        end
    end
    % Add legends
    if NUM_PROP + NUM_NODES > 2
        set_legends(ax, flip(legend_name), 'config', config)
    end
end
function gui_plot_custom_figures(obj)
    % Function that plot custom figures based on the parameters of the
    % GUI's custom figures tab
    
    % Get mixtures
    mixtures_nodes = obj.Tree_mixtures.CheckedNodes;
    % Get variables x
    x_field = get(obj.Tree_variable_x.CheckedNodes, 'Text');
    % Get variables y
    y_field = get(obj.Tree_variable_y.CheckedNodes, 'Text');
    if ~iscell(y_field)
        y_field = {y_field};
    end
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
    legend_name_1 = get(mixtures_nodes, 'Text')';
    legend_name_2 = y_field;
    % Clear axes
    cla(obj.UIAxes);
    % Initialize axes
    ax = set_figure(obj.UIAxes, config);
    % Plot
    k = NUM_NODES * NUM_PROP;
    for i = NUM_NODES:-1:1
        mixtures = obj.temp_results.(mixtures_nodes(i).Text);
        x_values = cell2vector(mixtures, x_field);
        for j = NUM_PROP:-1:1
            plot_figure(x_field, x_values, y_field{j}, mixtures, 'config', config, 'axes', obj.UIAxes, 'linestyle', LINE_STYLES{i}, 'color', COLOR_PALETTE(j, :));
            legend_name{k} = [legend_name_1{i}, ' - ', legend_name_2{j}];
            k = k - 1;
        end
    end
    % Add legends
    set_legends(ax, legend_name, config)
end
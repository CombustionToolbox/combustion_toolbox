function gui_plot_custom_figures(obj)
    % Function that plot custom figures based on the parameters of the
    % GUI's custom figures tab
    
    % Get mixtures
    mixtures_nodes = obj.Tree_mixtures.CheckedNodes;
    % Get variables x
    x_field = get(obj.Tree_variable_x.CheckedNodes, 'Text');
    % Get variables y
    y_field = get(obj.Tree_variable_y.CheckedNodes, 'Text');
    % Plot settings
    Misc = Miscellaneous();
    config = Misc.config;
    config.labelx = interpreter_label(x_field, config.label_type); % Set x label
    config.labely = interpreter_label(y_field, config.label_type); % Set y label
    % Clear axes
    cla(obj.UIAxes);
    % Initialize axes
    set_figure(obj.UIAxes, config);
    % Plot
    if ~iscell(y_field)
        y_field = {y_field};
    end
        
    for i = length(mixtures_nodes):-1:1
        mixtures = obj.temp_results.(mixtures_nodes(i).Text);
        x_values = cell2vector(mixtures, x_field);
        for j = length(y_field):-1:1
            plot_figure(x_field, x_values, y_field{j}, mixtures, 'config', config, 'axes', obj.UIAxes);
        end
    end
end
function gui_plot_custom_figures(obj)
    % Function that plot custom figures based on the parameters of the
    % GUI's custom figures tab
    
    % Get mixtures
    mixtures_nodes = obj.Tree_mixtures.CheckedNodes;
    % Get variables x
    x_name = get(obj.Tree_variable_x.CheckedNodes, 'Text');
    % Get variables y
    y_name = get(obj.Tree_variable_y.CheckedNodes, 'Text');
    % Plot settings
    obj.Misc.config.title = 'Figure';
    obj.Misc.config.labelx = interpret_label(x_name);
    obj.Misc.config.labely = interpret_label(y_name);
    obj.Misc.config.fontsize = 16;
    % Clear axes
    cla(obj.UIAxes);
    % Initialize axes
    set_figure(obj.UIAxes, obj.Misc.config);
    % Plot
    if ~iscell(y_name)
        y_name = {y_name};
    end
        
    for i = length(mixtures_nodes):-1:1
        mixtures = obj.temp_results.(mixtures_nodes(i).Text);
        x_values = cell2vector(mixtures, x_name);
        for j = length(y_name):-1:1
            plot_figure(x_values, mixtures, x_name, y_name{j}, obj.Misc.config, obj.PD.CompleteOrIncomplete, obj.UIAxes);
        end
    end
end
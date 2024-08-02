function gui_plot_custom_figures(app)
    % Function that plot custom figures based on the parameters of the
    % GUI's custom figures tab
    
    % Import packages
    import combustiontoolbox.utils.display.*
    import combustiontoolbox.utils.extensions.brewermap
    
    % Definitions
    basisName = 'mi';

    % Clear axes
    cla(app.UIAxes);

    % Get mixtures
    mixtures_nodes = app.Tree_mixtures.CheckedNodes;
    
    % Check name nodes
    nameNodes = {mixtures_nodes.Text};
    FLAG_REMOVE = strcmp(nameNodes, 'Mixtures');
    mixtures_nodes(FLAG_REMOVE) = [];

    % Get x and y variables
    x_field = get(app.Tree_variable_x.CheckedNodes, 'Text');
    y_field = get(app.Tree_variable_y.CheckedNodes, 'Text');

    % Check non empty values
    if isempty(x_field) || isempty(y_field), return; end

    % Check that variables y is a cell
    if ~iscell(y_field), y_field = {y_field}; end

    % Plot settings
    config = app.plotConfig();
    config.labelx = interpreterLabel(x_field, config.label_type); % Set x label
    config.labely = interpreterLabel(y_field, config.label_type); % Set y label
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
        ax = setFigure(app.UIAxes, config);
    end

    % Plot
    k = NUM_NODES * NUM_PROP;
    N = length(app.temp_results);

    for i = NUM_NODES:-1:1

        for j = N:-1:1
            mixtures{:, j} = app.temp_results(j).(mixtures_nodes(i).Text);
        end

        x_values = cell2vector(mixtures, x_field);
        basis = cell2vector(mixtures, basisName);

        x_values = checkUnits(x_field, x_values, basis);
        
        if any(strcmpi(y_field{1}, {'Xi', 'Yi', 'Ni', 'moles', 'molarFraction', 'massFraction'}))
            y_values = [mixtures{:}];
            plotComposition(y_values(1), y_values, x_field, y_field{1}, 'y_var', y_values, 'axes', app.UIAxes)
            return
        end

        for j = NUM_PROP:-1:1
            y_values = cell2vector(mixtures, y_field{j});
            y_values = checkUnits(y_field{j}, y_values, basis);
            plotFigure(x_field, x_values, y_field{j}, y_values, 'config', config, 'axes', app.UIAxes, 'linestyle', LINE_STYLES{i}, 'color', COLOR_PALETTE(j, :));
            legend_name{k} = [legend_name_1{i}, ' - ', legend_name_2{j}];
            k = k - 1;
        end
    end
    
    % Add legends and change ylabel
    if NUM_PROP + NUM_NODES > 2
        app.UIAxes.YLabel.String = 'Multiple variables';
        setLegends(app.UIAxes, flip(legend_name), 'config', config)
        app.UIAxes.Legend.Visible = 'on';
    else
        if isempty(app.UIAxes.Legend) 
            return
        end

        app.UIAxes.Legend.Visible = 'off';
    end

end

% SUB-PASS FUNCTIONS
function value = checkUnits(fieldname, value, valueBasis)
    % Change units if required
    switch lower(fieldname)
        case {'cp', 'cv', 'hf', 'ef', 'h', 'e', 'g', 's', 'det', 'dht', 'ds', 's0'}
            % Change units to [kJ ...]
            % value = value * 1e-3; 
            % Property has to be divided by the basis (kg or mol)
            value = value ./ valueBasis;
        case {'mw', 'w'}
            % Change units to [g ...]
            value = value * 1e3; 
    end

end
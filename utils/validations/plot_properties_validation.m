function mainfigure = plot_properties_validation(results1, results2, varsname_x, varsname_y, type, varargin)
    % Plot properties varname_y vs varname_x from CT (results1) against
    % results obtained from other code (results2). The properties to plot
    % are specified by varsname_x and varsname_y, which are cell arrays
    % of strings.
    
    % Default values
    nfrec = 1;
    config = results1.Misc.config;
    
    % Check plot type (with or without validation)
    FLAG_PLOT_VALIDATION = ~isempty(results2);

    % Create main figure
    mainfigure = figure;
    % set(mainfigure, 'position', config.position);
    set(mainfigure, 'units', 'normalized', 'innerposition', config.innerposition, 'outerposition', config.outerposition);
    
    tiledlayout(mainfigure, 'flow');

    % Loop
    for i = 1:length(varsname_x)
        varname_x = varsname_x{i};
        varname_y = varsname_y{i};

        dataname_x = get_dataname(varname_x, type);
        dataname_y = get_dataname(varname_y, type);
        results1.(varname_x) = cell2vector(select_data(results1, dataname_x), varname_x);
        results1.(varname_y) = cell2vector(select_data(results1, dataname_y), varname_y);

        nexttile;
        ax = gca;
        set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize - 2, 'BoxStyle', 'full')
        set(ax, 'TickLabelInterpreter', 'latex', 'Layer', 'top')
        grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on');
        xlabel(ax, interpreter_label(varname_x, config.label_type), 'FontSize', config.fontsize, 'interpreter', 'latex');
        ylabel(ax, interpreter_label(varname_y, config.label_type), 'FontSize', config.fontsize, 'interpreter', 'latex');
        xlim(ax, [min(results1.(varname_x)), max(results1.(varname_x))])

        NUM_COLORS = 1;
        LINE_STYLES = {'-', '--', ':', '-.'};
        SYMBOL_STYLES = {'d', 'o', 's', '<'};

        if NUM_COLORS > 1
            colorbw = brewermap(NUM_COLORS, config.colorpalette);
        else
            colorbw = config.colorlines(3, :);
        end

        k = 1;
        z = 1;

        plot(ax, results1.(varname_x), results1.(varname_y), 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'LineStyle', LINE_STYLES{z});

        if FLAG_PLOT_VALIDATION
            plot(ax, results2.(varname_x)(1:nfrec:end), results2.(varname_y)(1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'white');
        end

    end

end

% SUB-PASS FUNCTIONS
function dataname = get_dataname(var, type)

    if strcmpi(var, 'phi')
        dataname = 'PD.phi.value';
    elseif strcmpi(type, 'mix1')
        dataname = 'PS.strR';
    elseif strcmpi(type, 'mix2')
        dataname = 'PS.strP';
    elseif strcmpi(type, 'mix2_c')
        dataname = 'PS.mix2_c';
    elseif strcmpi(type, 'mix3')
        dataname = 'PS.mix3';
    end

end

function dataselected = select_data(self, dataname)
    index = strfind(dataname, '.');
    index = [index, length(dataname) + 1];
    N = length(index);
    dataselected = self;
    pos1 = 1;

    for i = 1:N
        pos2 = index(i) - 1;
        varname = dataname(pos1:pos2);
        dataselected = dataselected.(varname);
        pos1 = index(i) + 1;
    end

end
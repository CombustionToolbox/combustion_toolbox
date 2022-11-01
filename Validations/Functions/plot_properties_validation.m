function mainfigure = plot_properties_validation(results1, results2, varsname_x, varsname_y, type, varargin)
    % Plot properties varname_y vs varname_x from CT (results1) against results obtained from other code (results2)

    % Default values
    nfrec = 1;
    config = results1.Misc.config;
    % Create main figure
    mainfigure = figure;
    set(mainfigure, 'position', config.position);
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
        plot(ax, results2.(varname_x)(1:nfrec:end), results2.(varname_y)(1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'white');

        %         legendname = {'CT', 'CEA-NASA'};
        %         legend(ax, legendname, 'FontSize', config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
        %         title(create_title(results1), 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
    end

end

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

function titlename = create_title(self)
    titlename = [self.PD.ProblemType, ': ', cat_moles_species(self.PD.N_Fuel, self.PD.S_Fuel)];

    if ~isempty(self.PD.S_Oxidizer) || ~isempty(self.PD.S_Inert)

        if ~isempty(self.PD.S_Fuel)
            titlename = [titlename, ' + '];
        end

        titlename = [titlename, '$\frac{', sprintf('%.3g', self.PD.phi_t), '}{\phi}$'];

        if ~isempty(self.PD.S_Oxidizer) && ~isempty(self.PD.S_Inert)
            titlename = [titlename, '('];
        end

        if ~isempty(self.PD.S_Oxidizer)
            titlename = [titlename, cat_moles_species(1, self.PD.S_Oxidizer)];
        end

        if ~isempty(self.PD.S_Inert)
            titlename = [titlename, ' + ', cat_moles_species(self.PD.proportion_inerts_O2, self.PD.S_Inert)];
        end

        if ~isempty(self.PD.S_Oxidizer) && ~isempty(self.PD.S_Inert)
            titlename = [titlename, ')'];
        end

    end

end

function cat_text = cat_moles_species(moles, species)
    N = length(species);
    cat_text = [];

    if N
        cat_text = cat_mol_species(moles(1), species{1});

        for i = 2:N
            cat_text = [cat_text, ' + ', cat_mol_species(moles(i), species{i})];
        end

    end

end

function cat_text = cat_mol_species(mol, species)

    if mol == 1
        value = [];
    else
        value = sprintf('%.3g', mol);
    end

    cat_text = [value, species2latex(species)];
end

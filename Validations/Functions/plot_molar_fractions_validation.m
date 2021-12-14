function f = plot_molar_fractions_validation(results1, results2, varname_x, varname_y, species, varargin)
    % Default values
    nfrec = 1;
    mintol_display = 1e-16;
    config = results1.Misc.config;

    dataname_x = get_dataname(varname_x);
    dataname_y = get_dataname(varname_y);
    results1.(varname_x) = cell2vector(select_data(results1, dataname_x), varname_x);
    results1.(varname_y) = cell2vector(select_data(results1, dataname_y), varname_y);
    index_species_CT = find_ind(results1.Misc.LS_original, species);
    
    f = figure;
    set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
    axes = gca;
    set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')
    grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
    xlabel(axes, 'Equivalence ratio, $\phi$','FontSize',config.fontsize,'interpreter','latex');
    ylabel(axes, 'Molar fraction, $X_i$','FontSize',config.fontsize,'interpreter','latex');
%     set(axes,'xscale','log')
    set(axes,'yscale','log')
    xlim(axes, [min(results1.(varname_x)), max(results1.(varname_x))])
    ylim(axes, [mintol_display, 1])

    NE = length(species);
    maxLdisplay = config.colorpaletteLenght;
    if NE > maxLdisplay
        NUM_COLORS = maxLdisplay;
    else
        NUM_COLORS = NE;
    end

    LINE_STYLES = {'-', '--', ':', '-.'};
    SYMBOL_STYLES = {'d', 'o', 's', '<'};
    NUM_STYLES  = length(LINE_STYLES);
    colorbw = brewermap(NUM_COLORS, config.colorpalette);
    
    k = 1;
    z = 1;
    for i=1:length(species)
        plot(axes, results1.(varname_x), results1.(varname_y)(index_species_CT(i), :), 'LineWidth', config.linewidth, 'color', colorbw(k,:), 'LineStyle', LINE_STYLES{z});
        k = k + 1;
        if k == maxLdisplay
            k = 1;
            z = z + 1;
            if z > NUM_STYLES
                z = 1;
            end
        end
    end

    k = 1;
    z = 1;
    for i=1:length(species)
        plot(axes, results2.(varname_x)(1:nfrec:end), results2.(varname_y)(i, 1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k,:), 'MarkerFaceColor', 'white');
        k = k + 1;
        if k == maxLdisplay
            k = 1;
            z = z + 1;
            if z > NUM_STYLES
                z = 1;
            end
        end
    end
    
    k = 1;
    z = 1;
    h = zeros(1, length(species));
    for i=1:length(species)
        h(i) = plot(axes, NaN, NaN, strcat(LINE_STYLES{z}, SYMBOL_STYLES{z}), 'LineWidth', config.linewidth, 'color', colorbw(k,:));
        k = k + 1;
        if k == maxLdisplay
            k = 1;
            z = z + 1;
            if z > NUM_STYLES
                z = 1;
            end
        end
    end

    legendname = species;
    legend(h, legendname, 'FontSize', config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
    title(create_title(results1), 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
end

function dataname = get_dataname(var)
    if strcmpi(var, 'phi')
        dataname = 'PD.phi.value';
    else
        dataname = 'PS.strP';
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
    titlename = strcat(self.PD.ProblemType, ': ', cat_moles_species(self.PD.N_Fuel, self.PD.S_Fuel));
    if ~isempty(self.PD.S_Oxidizer) || ~isempty(self.PD.S_Inert)
        if ~isempty(self.PD.S_Fuel)
            titlename = strcat(titlename, ' + ');
        end
        titlename = strcat(titlename, strcat('$\frac{', sprintf('%.3g', self.PD.phi_t), '}{\phi}$'));
        if ~isempty(self.PD.S_Oxidizer) && ~isempty(self.PD.S_Inert)
            titlename = strcat(titlename, '(');
        end
        if ~isempty(self.PD.S_Oxidizer)
            titlename = strcat(titlename, cat_moles_species(1, self.PD.S_Oxidizer));
        end
        if ~isempty(self.PD.S_Inert)
            titlename = strcat(titlename, ' + ', cat_moles_species(self.PD.proportion_inerts_O2, self.PD.S_Inert));
        end
        if ~isempty(self.PD.S_Oxidizer) && ~isempty(self.PD.S_Inert)
            titlename = strcat(titlename, ')');
        end
    end
end

function cat_text = cat_moles_species(moles, species)
    N = length(species);
    cat_text = [];
    if N
        cat_text = cat_mol_species(moles(1), species{1});
        for i=2:N
            cat_text = strcat(cat_text, ' + ', cat_mol_species(moles(i), species{i}));
        end
    end
end

function cat_text = cat_mol_species(mol, species)
    if mol == 1
        value = [];
    else
        value = sprintf('%.3g', mol);
    end
    cat_text = strcat(value, convert_species_latex(species));
end

function speciesLatex = convert_species_latex(species)
    index = regexp(species, '[0-9]');
    N = length(index);
    if ~N
        speciesLatex = species;
        return;
    end
    speciesLatex = [];
    pos1 = 1;
    for i=1:N
        pos2 = index(i) - 1;
        speciesLatex = strcat(speciesLatex, species(pos1:pos2), '$_', species(pos2 + 1), '$');
        pos1 = pos2 + 2;
    end
    if pos1 < length(species)
        if isletter(species(pos1))
            aux = 0;
        else
            aux = 1;
        end
        speciesLatex = [speciesLatex, species(pos1+aux:end)];
    end
end
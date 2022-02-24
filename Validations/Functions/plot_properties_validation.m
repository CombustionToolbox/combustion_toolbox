function mainfigure = plot_properties_validation(results1, results2, varsname_x, varsname_y, varargin)
    % Plot properties varname_y vs varname_x from CT (results1) against CEA (results2)

    % Default values
    nfrec = 1;
    config = results1.Misc.config;
    % Create main figure
    mainfigure = figure;
    set(mainfigure, 'position', [1921 -471 1080 1795]);
    tiledlayout(mainfigure, 'flow');
    % Loop
    for i = 1:length(varsname_x)
        varname_x = varsname_x{i};
        varname_y = varsname_y{i};

        dataname_x = get_dataname(varname_x);
        dataname_y = get_dataname(varname_y);
        results1.(varname_x) = cell2vector(select_data(results1, dataname_x), varname_x);
        results1.(varname_y) = cell2vector(select_data(results1, dataname_y), varname_y);
    
        nexttile;
        axes = gca;
        set(axes, 'LineWidth', config.linewidth, 'FontSize', config.fontsize-2, 'BoxStyle', 'full')
        grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
        xlabel(axes, interpret_label(varname_x), 'FontSize', config.fontsize, 'interpreter', 'latex');
        ylabel(axes, interpret_label(varname_y), 'FontSize', config.fontsize, 'interpreter', 'latex');
        xlim(axes, [min(results1.(varname_x)), max(results1.(varname_x))])
    
        NUM_COLORS = 3;
        LINE_STYLES = {'-', '--', ':', '-.'};
        SYMBOL_STYLES = {'d', 'o', 's', '<'};
        colorbw = brewermap(NUM_COLORS, config.colorpalette);
        
        k = 1;
        z = 1;

        plot(axes, results1.(varname_x), results1.(varname_y), 'LineWidth', config.linewidth, 'color', colorbw(k,:), 'LineStyle', LINE_STYLES{z});
        plot(axes, results2.(varname_x)(1:nfrec:end), results2.(varname_y)(1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k,:), 'MarkerFaceColor', 'white');
        
%         legendname = {'CT', 'CEA-NASA'};
%         legend(axes, legendname, 'FontSize', config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
%         title(create_title(results1), 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
    end
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
    cat_text = strcat(value, species2latex(species));
end

function value = interpret_label(property)
    % Switch label from input property
    switch lower(property)
        case 'phi'
            value = 'Equivalance ratio';
        case 'rho'
            value = 'Density $[kg/m^3]$';
        case 't'
            value = 'Temperature $[K]$';
        case 'p'
            value = 'Pressure $[bar]$';
        case 'h'
            value = 'Enthalpy $[kJ/kg]$';
        case 'e'
            value = 'Internal energy $[kJ/kg]$';
        case 'g'
            value = 'Gibbs energy $[kJ/kg]$';
        case 's'
            value = 'Entropy $[kJ/kg-K]$';
        case 'cp'
            value = '$c_p [kJ/kg-K]$';
        case 'cv'
            value = '$c_v [kJ/kg-K]$';
        case 'gamma_s'
            value = 'Adiabatic index';
        otherwise
            value = '';
    end
end
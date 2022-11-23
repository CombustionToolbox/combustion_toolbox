function [ax, fig] = plot_molar_fractions(self, x_var, x_field, y_field, varargin)
    % Plot molar fractions againts any variable
    %
    % Args:
    %    self (struct): Data of the mixture, conditions, and databases
    %    x_var (cell): Properties of the mixture for all the cases
    %    x_field (char): Fieldname to plot on the x-axis
    %    y_field (char): Fieldname to plot on the y-axis
    %
    % Optional Args:
    %    * validation (struct): Struct that contains validations with (x_field, y_field)
    %    * nfrec (float): Frequency points to plot validations
    %    * mintol (float): Minimum limit i-axis with the composition of the mixture
    %    * config (struct): Struct with the configuration for plots
    %    * axis_x (char): Set x-axis limits
    %    * axis_y (char): Set y-axis limits
    %    * xscale (char): Set x-axis scale (linear or logarithmic)
    %    * yscale (char): Set y-axis scale (linear or logarithmic)
    %    * xdir (char): Set x-axis direction (normal or reverse)
    %    * ydir (char): Set y-axis direction (normal or reverse)
    %
    % Returns:
    %    Tuple containing
    %    * ax (axes): Axes object
    %    * fig (figure): Figure object

    % Default values
    ax = [];
    results2 = [];
    nfrec = 1;
    config = self.Misc.config;
    config.labelx = interpreter_label(x_field, config.label_type); % Set x label
    config.labely = interpreter_label(y_field, config.label_type); % Set y label
    mintol_display = self.C.mintol_display;
    config.yscale = 'log';
    [species, LS] = get_display_species(self);
    y_var = x_var;
    % Unpack
    for i = 1:2:nargin - 5

        switch lower(varargin{i})
            case {'validation', 'results'}
                results2 = varargin{i + 1};
            case {'nfrec'}
                nfrec = varargin{i + 1};
            case {'mintol', 'mintol_display', 'toln'}
                mintol_display = varargin{i + 1};
            case {'ls', 'species', 'display_species', 'display species'}
                species = varargin{i + 1};
            case {'y', 'y_var', 'yvar', 'y var', 'y_data', 'ydata', 'y data'}
                y_var = varargin{i + 1};
            case {'config'}
                config = varargin{i + 1};
            case {'xscale'}
                config.xscale = varargin{i + 1};
            case {'yscale'}
                config.yscale = varargin{i + 1};
            case {'xdir'}
                config.xdir = varargin{i + 1};
            case {'ydir'}
                config.ydir = varargin{i + 1};
            case {'ax', 'axes'}
                ax = varargin{i + 1};
        end

    end

    % Read data
    FLAG_Y_AXIS = strcmpi(y_field, 'Xi');
    if contains(self.PD.ProblemType, 'polar', 'IgnoreCase', true)
        mix1.(x_field) = cell2vector(x_var.polar, x_field);
        mix1.(y_field) = cell2vector(y_var.polar, y_field);
    else
        mix1.(x_field) = cell2vector(x_var, x_field);
        mix1.(y_field) = cell2vector(y_var, y_field);
    end
    % Get index species
    index_species_CT = find_ind(LS, species);
    % Remove species that do not appear
    [species, index_species_CT] = clean_display_species(mix1.Xi, LS, index_species_CT);
    % Set figure
    if isempty(ax)
        [ax, ~, fig] = set_figure(config);
    else
        fig = [];
    end

    % Set axis limits
    if FLAG_Y_AXIS
        xlim(ax, [min(mix1.(x_field)), max(mix1.(x_field))])
        ylim(ax, [mintol_display, 1])
    else
        xlim(ax, [mintol_display, 1])
        ylim(ax, [min(mix1.(y_field)), max(mix1.(y_field))])
    end

    % Set default style
    NE = length(species);
    maxLdisplay = config.colorpaletteLenght;

    if NE > maxLdisplay
        NUM_COLORS = maxLdisplay;
    else
        NUM_COLORS = NE;
    end

    LINE_STYLES = {'-', '--', ':', '-.'};
    SYMBOL_STYLES = {'d', 'o', 's', '<'};
    NUM_STYLES = length(LINE_STYLES);
    colorbw = brewermap(NUM_COLORS, config.colorpalette);
    % Plot main results
    k = 1;
    z = 1;

    for i = 1:length(species)

        if FLAG_Y_AXIS
            plot(ax, mix1.(x_field), mix1.(y_field)(index_species_CT(i), :), 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'LineStyle', LINE_STYLES{z});
        else
            plot(ax, mix1.(x_field)(index_species_CT(i), :), mix1.(y_field), 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'LineStyle', LINE_STYLES{z});
        end

        k = k + 1;

        if k == maxLdisplay
            k = 1;
            z = z + 1;

            if z > NUM_STYLES
                z = 1;
            end

        end

    end

    % Plot validations
    if ~isempty(results2)

        if iscell(results2)
            mix2.(x_field) = cell2vector(results2, x_field);
            mix2.(y_field) = cell2vector(results2, y_field);
        else
            mix2 = results2;
        end

        k = 1;
        z = 1;

        for i = 1:length(species)

            if FLAG_Y_AXIS
                plot(ax, mix2.(x_field)(1:nfrec:end), mix2.(y_field)(i, 1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'white');
            else
                plot(ax, mix2.(x_field)(i, 1:nfrec:end), mix2.(y_field)(1:nfrec:end), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'white');
            end

            k = k + 1;

            if k == maxLdisplay
                k = 1;
                z = z + 1;

                if z > NUM_STYLES
                    z = 1;
                end

            end

        end

        % Plot symbols
        k = 1;
        z = 1;
        h = zeros(1, length(species));

        for i = 1:length(species)
            h(i) = plot(ax, NaN, NaN, [LINE_STYLES{z}, SYMBOL_STYLES{z}], 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'auto');
            k = k + 1;

            if k == maxLdisplay
                k = 1;
                z = z + 1;

                if z > NUM_STYLES
                    z = 1;
                end

            end

        end

    else
        h = ax;
    end

    % Set legend
    for i = length(species):-1:1
        legendname{i} = species2latex(species{i});
    end

    legend(h, legendname, 'FontSize', config.fontsize - 6, 'Location', 'northeastoutside', 'interpreter', 'latex');
    % Set title
    title(create_title(self), 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
end

% SUB-PASS FUNCTIONS
function [species, LS] = get_display_species(self)
    LS = self.Misc.LS_original;

    if isempty(self.Misc.display_species)
        species = LS;
    else
        species = self.Misc.display_species;
    end

end

function [species, index_pass] = clean_display_species(molar_fractions, species, index_species)
    % Remove species that do not appear

    % Checks
    index_all_pass = find(any(molar_fractions > 0, 2));
    FLAG_PASS = ismember(index_species, index_all_pass);
    % Update species
    index_pass = index_species(FLAG_PASS);
    species = species(index_pass);
end

function titlename = create_title(self)
    label_problemtype = strrep(self.PD.ProblemType, '_', ' ');
    titlename = [label_problemtype, ': ', cat_moles_species(self.PD.N_Fuel, self.PD.S_Fuel)];

    if ~isempty(self.PD.S_Oxidizer) || ~isempty(self.PD.S_Inert)

        if ~isempty(self.PD.S_Fuel)
            titlename = [titlename, ' + '];
        end

        titlename = [titlename, '$\ \frac{', sprintf('%.3g', self.PD.phi_t), '}{\phi}$'];

        if ~isempty(self.PD.S_Oxidizer)

            if length(self.PD.S_Oxidizer) > 1
                titlename = [titlename, '('];
            end

            ind = find_ind(self.PD.S_Oxidizer, 'O2');

            if ind
                self.PD.N_Oxidizer = self.PD.N_Oxidizer / self.PD.N_Oxidizer(ind);
            end

            titlename = [titlename, cat_moles_species(self.PD.N_Oxidizer, self.PD.S_Oxidizer)];
        end

        if ~isempty(self.PD.S_Inert) && ~isempty(self.PD.ratio_inerts_O2)
            titlename = [titlename, ' + ', cat_moles_species(self.PD.N_Inert, self.PD.S_Inert)];
        end

        if ~isempty(self.PD.S_Oxidizer) && length(self.PD.S_Oxidizer) > 1
            titlename = [titlename, ')'];
        end

        if ~isempty(self.PD.S_Inert) && isempty(self.PD.ratio_inerts_O2)
            titlename = [titlename, ' + ', cat_moles_species(self.PD.N_Inert, self.PD.S_Inert)];
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

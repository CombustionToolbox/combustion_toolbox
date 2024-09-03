function [ax, fig] = plotComposition(obj, x_var, x_field, y_field, varargin)
    % Plot molar fractions againts any variable
    %
    % Args:
    %     obj (Mixture): Data of the mixture, conditions, and databases
    %     x_var (cell): Properties of the mixture for all the cases
    %     x_field (char): Field name for the x-axis data
    %     y_field (char): Field name for the y-axis data
    %
    % Optional Name-Value Pair Args:
    %     * validation (struct): Struct that contains validations with (x_field, y_field)
    %     * nfrec (float): Frequency points to plot validations
    %     * mintol (float): Minimum limit i-axis with the composition of the mixture
    %     * displaySpecies (cell): List of species to plot
    %     * y_var (cell): Get y-axis data from a different mixture
    %     * config (struct): Struct with the configuration for plots
    %     * axis_x (char): Set x-axis limits
    %     * axis_y (char): Set y-axis limits
    %     * xscale (char): Set x-axis scale (linear or log)
    %     * yscale (char): Set y-axis scale (linear or log)
    %     * xdir (char): Set x-axis direction (normal or reverse)
    %     * ydir (char): Set y-axis direction (normal or reverse)
    %     * ax (object): Handle of the axes to plot on
    %
    % Returns:
    %     Tuple containing
    %
    %     * ax (object): Handle of the axes
    %     * fig (object): Handle of the figure
    %
    % Examples:
    %     * [ax, fig] = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi')
    %     * [ax, fig] = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'y_var', mixArray2)
    %     * [ax, fig] = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'y_var', mixArray2, 'validation', resultsCEA)
    %     * [ax, fig] = plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'y_var', mixArray2, 'validation', resultsCEA, 'displaySpecies', displaySpecies)
    
    % Import packages
    import combustiontoolbox.utils.findIndex
    import combustiontoolbox.utils.display.*
    
    % Initialization
    config = PlotConfig();
    mintolDisplay = config.mintolDisplay;
    displaySpecies = config.displaySpecies;

    % Default values
    ax = [];
    results2 = [];
    nfrec = 1;
    FLAG_PLOT_VALIDATION = false;
    % config = self.Misc.config;
    config.labelx = interpreterLabel(x_field, config.label_type);
    config.labely = interpreterLabel(y_field, config.label_type);
    % mintol_display = self.C.mintol_display;
    config.yscale = 'log';
    y_var = x_var;

    % Unpack
    for i = 1:2:nargin - 5

        switch lower(varargin{i})
            case {'validation', 'results'}
                results2 = varargin{i + 1};
                FLAG_PLOT_VALIDATION = true;
            case {'nfrec'}
                nfrec = varargin{i + 1};
            case {'mintol', 'mintol_display', 'toln'}
                mintolDisplay = varargin{i + 1};
            case {'ls', 'species', 'displayspecies', 'display species', 'display_species'}
                displaySpecies = varargin{i + 1};
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
            case {'title'}
                config.title = varargin{i + 1};
            case {'ax', 'axes'}
                ax = varargin{i + 1};
        end

    end
    
    % Set species to displaySpecies
    [species, listSpecies] = get_displaySpecies(obj, displaySpecies);

    % Read data
    FLAG_Y_AXIS = strcmpi(y_field, 'Xi');
    % if contains(self.PD.ProblemType, 'polar', 'IgnoreCase', true)
    %     mix1.(x_field) = cell2vector(x_var.polar, x_field);
    %     mix1.(y_field) = cell2vector(y_var.polar, y_field);
    % elseif ~(iscell(x_var) && iscell(y_var))
    %     mix1.(x_field) = x_var;
    %     mix1.(y_field) = y_var;
    % else
    %     mix1.(x_field) = cell2vector(x_var, x_field);
    %     mix1.(y_field) = cell2vector(y_var, y_field);
    % end
    if isfloat(x_var) && isfloat(y_var)
        mix1.(x_field) = x_var;
        mix1.(y_field) = y_var;
    elseif ~(iscell(x_var) && iscell(y_var))
        mix1.(x_field) = [x_var.(x_field)];
        mix1.(y_field) = [y_var.(y_field)];
    else
        mix1.(x_field) = cell2vector(x_var, x_field);
        mix1.(y_field) = cell2vector(y_var, y_field);
    end
    
    % Check lenghts
    if length(mix1.(x_field)) < 2
        return
    end

    % Get index species
    index_species_CT = findIndex(listSpecies, species);

    % Remove species that do not appear
    [species, index_species_CT] = clean_displaySpecies(mix1.Xi, listSpecies, index_species_CT, mintolDisplay, FLAG_PLOT_VALIDATION);

    % Set figure
    if isempty(ax)
        [ax, ~, fig] = setFigure(config);
    else
        fig = [];
    end

    % Set axis limits
    if FLAG_Y_AXIS
        xlim(ax, [min(mix1.(x_field)), max(mix1.(x_field))])
        ylim(ax, [mintolDisplay, 1])
    else
        xlim(ax, [mintolDisplay, 1])
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
    colorbw = combustiontoolbox.utils.extensions.brewermap(NUM_COLORS, config.colorpalette);

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
    if FLAG_PLOT_VALIDATION

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

    legend(h, legendname, 'FontSize', config.fontsize - 2, 'Location', 'northeastoutside', 'interpreter', 'latex');
    
    % Set title
    if isempty(config.title)
        config.title = getTitle(obj);
    else
        config.title = sprintf('%s - %s', config.title, getTitle(obj));
    end
    
    title(config.title, 'Interpreter', 'latex', 'FontSize', config.fontsize + 4);
end

% SUB-PASS FUNCTIONS
function [species, listSpecies] = get_displaySpecies(obj, displaySpecies)
    listSpecies = obj.chemicalSystem.listSpecies;

    if isempty(displaySpecies)
        species = obj.chemicalSystem.listProducts;
    else
        species = displaySpecies;
    end

end

function [species, index_pass] = clean_displaySpecies(molar_fractions, species, index_species, mintolDisplay, FLAG_PLOT_VALIDATION)
    % Remove species that do not appear

    % Checks
    if FLAG_PLOT_VALIDATION
        mintolDisplay = 0;
    end

    index_all_pass = find(any(molar_fractions > mintolDisplay, 2));
    FLAG_PASS = ismember(index_species, index_all_pass);
    % Update species
    index_pass = index_species(FLAG_PASS);
    species = species(index_pass);
end
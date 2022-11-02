function ax = plot_thermo_validation(species, property, DB, varargin)
    % Validation custom thermodynamic polynomials with NASA's 9 polynomials
    %
    % Args:
    %    species (cell): List of species
    %    property (str): Name of the thermodynamic property to check
    %    DB (struct):    Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Optional Args:
    %    nfrec(float):   Points frequency for NASA values
    %    range(float):   Temperature range [K]
    %
    % Returns:
    %    ax (axes): Axes of the plotted figure

    % Check input species
    if ischar(species)
        species = {species};
    end

    % Definitions
    NS = length(species);
    nfrec = 1;
    x_labelname = 'Temperature T [K]';
    FLAG_RANGE = true;
    % Unpack inputs
    nargin0 = 3;

    for i = nargin - nargin0:-2:1

        switch lower(varargin{i - 1})
            case 'nfrec'
                nfrec = varargin{i};
            case 'range'
                range = varargin{i};
                range = repmat(range, NS, 1)';
                FLAG_RANGE = false;
        end

    end

    % Set input values
    [funname_NASA, funname_CT, y_labelname] = set_inputs_thermo_validations(property);
    % Compute values
    for i = NS:-1:1
        % NASA's polynomials
        if FLAG_RANGE
            range(:, i) = DB.(species{i}).T;
        end

        fun_NASA = @(T) funname_NASA(species{i}, T, DB);
        result_NASA(:, i) = feval(fun_NASA, range(:, i));
        % Combustion Toolbox
        fun_CT = @(T) funname_CT(species{i}, T, DB);
        result_CT(:, i) = feval(fun_CT, range(:, i));
    end

    % Plot results
    ax = plot_results(species, NS, range, result_CT, result_NASA, x_labelname, y_labelname, nfrec);
end

% SUB-PASS FUNCTIONS
function ax = plot_results(species, NS, range, result_CT, result_NASA, x_labelname, y_labelname, nfrec)
    % Plot results

    [ax, config] = set_figure();
    set(ax, 'XScale', 'log')

    maxLdisplay = config.colorpaletteLenght;

    if NS > maxLdisplay
        NUM_COLORS = maxLdisplay;
    else
        NUM_COLORS = NS;
    end

    LINE_STYLES = {'-', '--', ':', '-.'};
    SYMBOL_STYLES = {'d', 'o', 's', '<'};
    NUM_STYLES = length(LINE_STYLES);
    colorbw = brewermap(NUM_COLORS, 'Dark2');

    k = 1;
    z = 1;

    for i = 1:NS
        plot(ax, range(:, i), result_CT(:, i), 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'LineStyle', LINE_STYLES{z});
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

    for i = 1:NS
        plot(ax, range(1:nfrec:end, i), result_NASA(1:nfrec:end, i), SYMBOL_STYLES{z}, 'LineWidth', config.linewidth, 'color', colorbw(k, :), 'MarkerFaceColor', 'white');
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
    h = zeros(1, NS);

    for i = 1:NS
        h(i) = plot(ax, NaN, NaN, strcat(LINE_STYLES{z}, SYMBOL_STYLES{z}), 'LineWidth', config.linewidth, 'color', colorbw(k, :));
        k = k + 1;

        if k == maxLdisplay
            k = 1;
            z = z + 1;

            if z > NUM_STYLES
                z = 1;
            end

        end

    end

    for i = NS:-1:1
        legendname{i} = species2latex(species{i});
    end

    legend(h, legendname, 'FontSize', config.fontsize - 2, 'Location', 'northeastoutside', 'interpreter', 'latex');

    xlabel(ax, x_labelname)
    ylabel(ax, y_labelname)
end

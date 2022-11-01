function ax = displaysweepresults(self, mix, xvar, varargin)
    % Plot a given variable against the molar fractions of the mixture
    %
    % Args:
    %    self (struct): Data of the mixture, conditions, and databases
    %    mix (struct): Properties of the mixture
    %    xvar (float): Vector with the x data
    %
    % Returns:
    %    ax (axes): Axes object

    % Default values
    config = self.Misc.config;
    display_species = self.Misc.display_species;
    species = self.S.LS;
    mintol = self.C.mintol_display;
    xscale = 'linear';
    yscale = 'log';
    xdir = 'normal';
    ydir = 'normal';
    % Unpack
    for i = 1:2:nargin - 3

        switch lower(varargin{i})
            case {'config'}
                config = varargin{i + 1};
            case {'xscale'}
                xscale = varargin{i + 1};
            case {'yscale'}
                yscale = varargin{i + 1};
            case {'xdir'}
                xdir = varargin{i + 1};
            case {'ydir'}
                ydir = varargin{i + 1};
        end

    end

    % Function
    if length(xvar) > 1
        ax = set_figure(config);
        set(ax, 'LineWidth', config.linewidth, 'FontSize', config.fontsize - 2, 'BoxStyle', 'full')
        grid(ax, 'off'); box(ax, 'off'); hold(ax, 'on'); ax.Layer = 'Top';
        set(ax, 'xscale', xscale, 'yscale', yscale)
        set(ax, 'XDir', xdir, 'YDir', ydir)
        [yvar, all_ind, indlabel] = get_parameters(mix, species, display_species, mintol);
        plot_line(self, ax, xvar, yvar, all_ind, indlabel, species)
        set_limits_axes(ax, xvar, yvar, mintol);
        species = species(all_ind);

        for i = length(species):-1:1
            legendname{i} = species2latex(species{i});
        end

        set_legends(self, ax, legendname);
    end

end

% SUB-PASS FUNCTIONS
function plot_line(self, ax, xvar, yvar, indy, indlabel, label_name)
    NE = numel(indy);
    maxLdisplay = self.Misc.config.colorpaletteLenght;

    if NE > maxLdisplay
        NUM_COLORS = maxLdisplay;
    else
        NUM_COLORS = NE;
    end

    LINE_STYLES = {'-', '--', ':', '-.'};
    NUM_STYLES = length(LINE_STYLES);

    colorbw = brewermap(NUM_COLORS, self.Misc.config.colorpalette);

    k = 1;
    z = 1;

    for i = 1:numel(indy)
        dl = plot(ax, xvar, yvar(indy(i), :), 'LineWidth', self.Misc.config.linewidth, 'color', colorbw(k, :), 'LineStyle', LINE_STYLES{z});
        label_line(self, ax, dl, label_name{indlabel(i)}, i)
        k = k + 1;

        if k == maxLdisplay
            k = 1;
            z = z + 1;

            if z > NUM_STYLES
                z = 1;
            end

        end

    end

end

function label_line(self, ax, dl, label_name, i)

    if self.Misc.FLAG_LABELS

        if mod(i, nfrec) == 0
            loc_label = 'right';
        else
            loc_label = 'left';
        end

        label(dl, label_name, 'FontSize', self.Misc.config.fontsize, ...
            'location', loc_label, 'ax', ax);
    end

end

function set_limits_axes(ax, xvar, yvar, mintol)
    set(ax, 'yscale', 'log')
    xlim(ax, [min(xvar), max(xvar)])
    ymin = 10^floor(log(abs(min(yvar(yvar > 0)))) / log(10));

    if ymin > mintol
        ylim(ax, [ymin, 1])
    else
        ylim(ax, [mintol, 1])
    end

end

function set_legends(self, ax, legends_name)
    legend(ax, legends_name, 'FontSize', self.Misc.config.fontsize - 6, 'Location', 'northeastoutside', 'Interpreter', 'latex');
end

function [yvar, all_ind, indlabel] = get_parameters(mix, species, display_species, mintol)

    if isfield(mix, 'polar')
        [yvar, all_ind, indlabel] = get_parameters_polar(mix.polar, species, display_species, mintol);
    else
        [yvar, all_ind, indlabel] = get_parameters_default(mix, species, display_species, mintol);
    end

end

function [yvar, all_ind, indlabel] = get_parameters_default(mix, species, display_species, mintol)
    Nstruct = length(mix);
    yvar = zeros(length(mix{1}.Xi), Nstruct);
    all_ind = [];

    if isempty(display_species)

        for i = Nstruct:-1:1
            yvar(:, i) = mix{i}.Xi;
            j = mix{i}.Xi > mintol;
            ind = find(j > 0);
            all_ind = [all_ind; ind];
        end

        all_ind = unique(all_ind);
        indlabel = all_ind;
    else

        for i = Nstruct:-1:1
            yvar(:, i) = mix{i}.Xi;
        end

        for i = 1:length(species)

            if any(strcmpi(species{i}, display_species))
                all_ind = [all_ind, i];
            end

        end

        indlabel = any(yvar(all_ind, :)' > mintol);
        all_ind = all_ind(indlabel);
    end

end

function [yvar, all_ind, indlabel] = get_parameters_polar(mix, species, display_species, mintol)
    Nstruct = length(mix.theta);
    yvar = zeros(size(mix.Xi));
    all_ind = [];

    if isempty(display_species)

        for i = Nstruct:-1:1
            yvar(:, i) = mix.Xi(:, i);
            j = mix.Xi(:, i) > mintol;
            ind = find(j > 0);
            all_ind = [all_ind; ind];
        end

        all_ind = unique(all_ind);
        indlabel = all_ind;
    else

        for i = Nstruct:-1:1
            yvar(:, i) = mix.Xi(:, i);
        end

        for i = 1:length(species)

            if any(strcmpi(species{i}, display_species))
                all_ind = [all_ind, i];
            end

        end

        indlabel = any(yvar(all_ind, :)' > mintol);
        all_ind = all_ind(indlabel);
    end

end

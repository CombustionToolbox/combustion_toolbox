function post_results(self)
    % Postprocess all the results with predefined plots
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases

    % Export results
    export_results(self);
    
    % Get problem type
    ProblemType = self.PD.ProblemType;
    % Get value attribute parametric study
    range = self.PD.range;
    % Get FLAG_RANGE == parametric study
    FLAG_PLOT_RANGE = numel(range) > 1 && all(range(2:end) ~= range(1));
    % Check FLAG_RANGE
    if ~FLAG_PLOT_RANGE && ~contains(ProblemType, 'POLAR', 'IgnoreCase', true)
        return
    end

    % Initialization
    [self, mix2, range_name, report_type] = plot_initialize(self);
    
    % Plot products composition vs. attribute of the parametric study
    switch lower(ProblemType)
        case {'shock_polar', 'det_polar'}
            % Post-shock
            plot_composition_polar(self, mix2)

        case 'shock_polar_r'
            % Post-shock - incident
            plot_composition_polar(self, self.PS.str2, 'Incident')
            % Post-shock - reflected
            plot_composition_polar(self, mix2, 'Reflected')
        
        case {'shock_r', 'det_r', 'det_overdriven_r', 'det_underdriven_r', 'det_oblique'}
            % Post-shock - incident
            plot_composition(self, self.PS.str2, range_name, range, 'Incident')
            % Post-shock - reflected
            plot_composition(self, mix2, range_name, range, 'Reflected')

        otherwise
            plot_composition(self, mix2, range_name, range);
    end
    
    % Select type of report
    switch lower(report_type)
        case 'short'
            plot_properties_short(self, range_name, range, ProblemType);
        case 'complete'
            plot_properties_complete(self, range_name, range, ProblemType);
    end

end

% SUB-PASS FUNCTIONS
function [self, mix2, range_name, report_type] = plot_initialize(self)
    % Initialization routine

    % Get label attribute of the parametric study
    range_name = self.PD.range_name;
    % Remove suffix
    range_name(range_name == 'R' | range_name == 'P') = [];
    % Remove subscripts
    range_name(range_name == '1') = [];
    % Set title
    self.Misc.config.title = get_title(self);
    % Type of report
    report_type = self.Misc.report_type;

    % Get mixtures
    try
        mix2 = self.PS.strP;
    catch
        mix2 = self.PS.str2;
    end
    
    if isempty(mix2{end})
        mix2 = self.PS.mix3;
    end

end

function dataname = get_dataname(mix_label)
    % Get data label
    
    switch lower(mix_label)
        case {'mix1', 'strr'}
            dataname = 'PS.strR';
        case {'mix2', 'strp'}
            dataname = 'PS.strP';
        case {'mix2_c'}
            dataname = 'PS.mix2_c';
        case {'mix3'}
            dataname = 'PS.mix3';
        case {'mix4'}
            dataname = 'PS.mix4';
        otherwise
            dataname = sprintf('PS.%s', mix_label);
    end

end

function dataselected = select_data(self, dataname)
    % Get data selected

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

function plot_composition(self, mix, range_name, range, varargin)
    % Post-processing routine that plots the molar composition againts the
    % given vector (range)
    
    % Defaults
    prefix = [];
    composition_label = 'Xi';
    
    % Unpack
    if nargin > 4
        prefix = varargin{1};
        prefix = [prefix, ' | '];
    end

    

    % Get composition
    composition = cell2vector(mix, composition_label);

    % Plot composition
    ax = plot_molar_fractions(self, range, range_name, composition_label, 'yvar', composition);
    
    % Update title
    if isempty(prefix)
        return
    end

    self.Misc.config.title = [prefix, self.Misc.config.title];
    set_title(ax, self.Misc.config)
end

function plot_composition_polar(self, mix, varargin)
    % Post-processing routine that plots the molar composition againts the
    % given vector (range)
    
    % Defaults
    prefix = [];
    composition_label = 'Xi';
    
    % Unpack
    if nargin > 2
        prefix = varargin{1};
        prefix = [prefix, ' | '];
    end

    % Definitions
    M1 = cell2vector(self.PS.strR, 'u') ./ cell2vector(self.PS.strR, 'sound');
    title_0 = self.Misc.config.title;

    % Plot composition
    for i = 1:length(mix)
        self.Misc.config.title = [prefix, '$\mathcal{M}_1 ', sprintf('= %.2f$ | %s', M1(i)), title_0];
        ax = plot_molar_fractions(self, mix{i}, 'beta', composition_label);
        set_title(ax, self.Misc.config)
    end

end

function plot_properties_short(self, range_name, range, ProblemType)
    % Post-processing routine that plots a set of properties againts the
    % given vector (range)
    
    % Default
    FLAG_RETURN = false;
    FLAG_PLOT = true;
    title_plots = {'Products'};
    properties_basis = [];

    % Set properties to plot
    switch lower(ProblemType)
        case {'tp', 'tv'}
            return

        case {'hp', 'ev', 'sp', 'sv'}
            mix_labels = {'mix2'};
            properties = {'T'};

        case {'shock_i'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_PLOT = [false, true];
            title_plots = {'Post-shock'};
            properties = {'T', 'p'};

        case {'shock_r', 'det_r', 'det_overdriven_r', 'det_underdriven_r'}
            mix_labels = {'mix1', 'str2', 'mix2'};
            FLAG_PLOT = [false, true, true];
            title_plots = {'Incident', 'Reflected'};
            properties = {'T', 'p'};
        
        case {'shock_oblique'}
            if strcmpi(range_name, 'theta')
                mix_labels = {'mix2', 'str2'};
                title_plots = {'Strong shock', 'Weak shock'};
                FLAG_PLOT = true(1, 2);
            else
                mix_labels = {'mix2'};
                title_plots = {'Post-shock'};
            end

            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};

        case {'shock_polar', 'det_polar'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_RETURN = true;

        case {'shock_polar_r'}
            mix_labels = {'mix1', 'str2', 'str2_1', 'mix2'};
            title_plots = create_cell_ntimes('Products', 4);
            FLAG_PLOT = true(1, 4);
            FLAG_RETURN = true;
        
        case {'det', 'det_overdriven', 'det_underdriven'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_PLOT = [false, true];
            title_plots = {'Post-detonation'};
            properties = {'T', 'p'};
        
        case {'det_oblique'}
            mix_labels = {'mix2', 'str2'};
            title_plots = {'Strong detonation', 'Weak detonation'};
            FLAG_PLOT = true(1, 2);

            properties = {'T', 'p'};

        case {'rocket'}
            mix_labels = {'mix2'};
            properties = {'T', 'p', 'I_sp', 'I_vac'};
    end
    
    % Definitions
    Nmix = length(mix_labels);
    
    % Get mixtures
    for i = Nmix:-1:1
        dataname{i} = get_dataname(mix_labels{i});
        mix{:, i} = select_data(self, dataname{i});
    end
    
    % First plot non-common plots
    plot_additional(self, mix, ProblemType);
    
    % Check FLAG_RETURN
    if FLAG_RETURN
        return
    end
    
    % Remove mixtures that have not to be plotted
    mix(:, ~FLAG_PLOT) = [];

    % Miscellaneous
    title_0 = self.Misc.config.title;

    % Plot properties vs. range
    for i = 1:length(mix)
        self.Misc.config.title = sprintf('%s | %s', title_plots{i}, title_0);
        plot_figure(range_name, range, properties, mix{:, i}, 'config', self.Misc.config, 'basis', properties_basis);
    end

end

function plot_properties_complete(self, range_name, range, ProblemType)
    % Post-processing routine that plots a set of properties againts the
    % given vector (range)
    
    % Default
    FLAG_RETURN = false;
    FLAG_PLOT = true;
    title_plots = {'Products'};
    properties_basis = [];

    % Set properties to plot
    switch lower(ProblemType)
        case {'tp', 'tv', 'hp', 'ev', 'sp', 'sv'}
            mix_labels = {'mix2'};
            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};

        case {'shock_i'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_PLOT = [false, true];
            title_plots = {'Post-shock'};
            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};

        case {'shock_r', 'det_r', 'det_overdriven_r', 'det_underdriven_r'}
            mix_labels = {'mix1', 'str2', 'mix2'};
            FLAG_PLOT = [false, true, true];
            title_plots = {'Incident', 'Reflected'};
            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};
        
        case {'shock_oblique'}
            if strcmpi(range_name, 'theta')
                mix_labels = {'mix2', 'str2'};
                title_plots = {'Strong shock', 'Weak shock'};
                FLAG_PLOT = true(1, 2);
            else
                mix_labels = {'mix2'};
                title_plots = {'Post-shock'};
            end

            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};

        case {'shock_polar', 'det_polar'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_RETURN = true;

        case {'shock_polar_r'}
            mix_labels = {'mix1', 'str2', 'str2_1', 'mix2'};
            title_plots = create_cell_ntimes('Products', 4);
            FLAG_PLOT = true(1, 4);
            FLAG_RETURN = true;
        
        case {'det', 'det_overdriven', 'det_underdriven'}
            mix_labels = {'mix1', 'mix2'};
            FLAG_PLOT = [false, true];
            title_plots = {'Post-detonation'};
            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};
        
        case {'det_oblique'}
            mix_labels = {'mix2', 'str2'};
            title_plots = {'Strong detonation', 'Weak detonation'};
            FLAG_PLOT = true(1, 2);

            properties = {'T', 'p', 'rho', 'h', 'e', 'g', 'S', 'gamma_s'};
            properties_basis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};

        case {'rocket'}
            mix_labels = {'mix2'};
            properties = {'T', 'p', 'rho', 'h', 'gamma_s', 'cstar', 'I_sp', 'I_vac'};
            properties_basis = {[], [], [], 'mi', [], [], [], []};
    end
    
    % Definitions
    Nmix = length(mix_labels);
    
    % Get mixtures
    for i = Nmix:-1:1
        dataname{i} = get_dataname(mix_labels{i});
        mix{:, i} = select_data(self, dataname{i});
    end
    
    % First plot non-common plots
    plot_additional(self, mix, ProblemType);
    
    % Check FLAG_RETURN
    if FLAG_RETURN
        return
    end
    
    % Remove mixtures that have not to be plotted
    mix(:, ~FLAG_PLOT) = [];

    % Miscellaneous
    title_0 = self.Misc.config.title;

    % Plot properties vs. range
    for i = 1:length(mix)
        self.Misc.config.title = sprintf('%s | %s', title_plots{i}, title_0);
        plot_figure_set(range_name, range, properties, mix{:, i}, 'config', self.Misc.config, 'basis', properties_basis);
    end

end

function plot_additional(self, mix, ProblemType)
    % Additional plots

    switch lower(ProblemType)
        case {'shock_i', 'det', 'det_underdriven', 'det_overdriven'}
            % Plot Hugoniout curve
            plot_hugoniot(self, mix{:, 1}, mix{:, 2});
        
        case {'shock_r', 'det_r', 'det_underdriven_r', 'det_overdriven_r'}
            % Plot Hugoniout curve - incident
            ax = plot_hugoniot(self, mix{:, 1}, mix{:, 2});
            % Plot Hugoniout curve - reflected
            ax = plot_hugoniot(self, mix{:, 1}, mix{:, 3}, ax);
            % Set legends
            self.Misc.config.legend_name = {'incident', 'reflected'};
            set_legends(ax, self.Misc.config.legend_name, 'config', self.Misc.config);

        case {'shock_polar', 'det_polar'}
            % Shock/detonation polar diagrams - incident
            plot_shock_polar(self, mix{:, 1}, mix{:, 2});

        case {'shock_polar_r'}
            % Shock polar diagrams - incident
            plot_shock_polar(self, mix{:, 1}, mix{:, 2});
            % Shock polar diagrams - reflected
            plot_shock_polar(self, mix{:, 2}, mix{:, 4}, mix{:, 3}, mix{:, 1});

    end

end
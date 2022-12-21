function post_results(self)
    % Postprocess all the results with predefined plots
    %
    % Args:
    %    self (struct): Data of the mixture, conditions, and databases

    % Abbreviations ---------------------
    mix1 = self.PS.strR;
    try
        mix2 = self.PS.strP;
    catch
        mix2 = self.PS.str2;
    end
    phi = self.PD.phi.value;
    range = self.PD.range;
    range_name = self.PD.range_name;
    ProblemType = self.PD.ProblemType;
    % -----------------------------------
    titlename = strrep(ProblemType, '_', ' ');
    self.Misc.config.title = titlename;
    FLAG_PLOT_RANGE = numel(range) > 1 && all(range(2:end) ~= range(1));
    FLAG_PLOT_PHI = numel(phi) > 1 && all(phi(2:end) ~= phi(1));

    if ~isfield(self.Misc.FLAGS_PROP, 'TP')
        self.Misc.FLAGS_PROP.TP = false;
    end

    if isempty(mix2{end})
        mix2 = self.PS.mix3;
    end

    % PLOTS REACTANTS
    if (self.Misc.FLAGS_PROP.TR || self.Misc.FLAGS_PROP.pR) && FLAG_PLOT_RANGE
        plot_molar_fractions(self, mix2, self.PD.range_name, 'Xi');
    elseif isfield(self.Misc.FLAGS_PROP, 'pP')

        if self.Misc.FLAGS_PROP.pP && FLAG_PLOT_RANGE
            plot_molar_fractions(self, mix2, self.PD.range_name, 'Xi');
        end

    end

    %     OF = cell2vector(self.PS.strR, 'OF');

    % PLOTS PRODUCTS
    if FLAG_PLOT_PHI
        plot_molar_fractions(self, mix1, 'phi', 'Xi', 'ydata', mix2);

        if ~any(strcmp(ProblemType, {'TP', 'TV'}))
            plot_figure('phi', range, 'T', mix2, 'config', self.Misc.config);
        end

    end

    if self.Misc.FLAGS_PROP.TP && length(range) > 1
        plot_molar_fractions(self, mix2, 'T', 'Xi');

    elseif strcmp(ProblemType, {'SP'}) && length(range) > 1
        self.Misc.config.labelx = 'Pressure $p$ [bar]';
        self.Misc.config.labely = 'Temperature $T$ [K]';
        plot_figure('p', self.PD.pP.value, 'T', mix2, 'config', self.Misc.config);

    elseif strcmp(ProblemType, {'SV'}) && length(range) > 1
        self.Misc.config.labelx = 'Volume ratio $v_P/v_R$';
        self.Misc.config.labely = 'Temperature $T$ [K]';
        plot_figure('v_P/v_R', self.PD.vP_vR.value, 'T', mix2, 'config', self.Misc.config);
        plot_molar_fractions(self, mix2, 'T', 'Xi');
    elseif strcmp(ProblemType, {'SHOCK_POLAR'}) || strcmp(ProblemType, {'DET_POLAR'})
        % Shock polars
        plot_shock_polar(self, mix1, mix2);
        %     % Molar fractions [-] againts incident velocity [m/s]
        %     self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
        %     self.Misc.config.labely = 'Molar fraction $X_i$';
        %     for i = length(mix2):-1:1
        %         beta = cell2vector(mix2{i}.polar, 'beta');
        %         displaysweepresults(self, mix2{i}.polar, beta);
        %     end
    elseif strcmp(ProblemType, {'SHOCK_POLAR_R'})
        mix3 = mix2;
        mix2 = self.PS.str2;
        mix2_1 = self.PS.str2_1;
        % mix3_1 = self.PS.str3_1;
        % mix3_2 = self.PS.str3_2;
        % Shock polars incident
        plot_shock_polar(self, mix1, mix2);
        % Shock polars reflected
        plot_shock_polar(self, mix2, mix3, mix2_1, mix1);
        % molar fractions [-] againts wave angle [deg]
        self.Misc.config.labelx = 'Wave angle $\beta$ [deg]';
        self.Misc.config.labely = 'Molar fraction $X_i$';

        for i = length(mix2):-1:1
            titlename = 'Deflection angle $\theta';
            self.Misc.config.title = sprintf('%s = %.2f$ [deg]', titlename, mix2_1{i}.theta);
            plot_molar_fractions(self, mix3{i}, 'beta', 'Xi');
            %         xlim(ax, [min(mix2{i}.polar.beta), 90]);
        end

%         % Temperature [K] againts deflection angle [deg]
%         titlename = '\mathcal{M}_1';
%         self.Misc.config.title = sprintf('%s = %.2f', titlename, mix1{i}.u / mix1{i}.sound);
%         ax = plot_figure('theta', mix2_1{1}.polar, 'T', mix2_1{1}.polar, 'config', self.Misc.config);
%         plot_figure('theta', mix3_1, 'T', mix3_1, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
%         plot_figure('theta', mix3_2, 'T', mix3_2, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
% 
%         ax = plot_figure('theta', mix2_1, 'p', mix2_1, 'config', self.Misc.config);
%         plot_figure('theta', mix3_1, 'p', mix3_1, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
%         plot_figure('theta', mix3_2, 'p', mix3_2, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
    elseif strcmp(ProblemType, {'SHOCK_I'}) && length(range) > 1
        self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
        self.Misc.config.labely = 'Temperature $T$ [K]';
        plot_figure('u', self.PD.u1.value, 'T', mix2, 'config', self.Misc.config);
        plot_hugoniot(self, mix1, mix2);
        self.Misc.config.labelx = 'Temperature $T$ [K]';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        T = cell2vector(mix2, 'T');
        displaysweepresults(self, mix2, T);
        self.Misc.config.labelx = 'Pressure $p$ [bar]';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        p = cell2vector(mix2, 'p');
        displaysweepresults(self, mix2, p);
        self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        u = cell2vector(mix1, 'u');
        displaysweepresults(self, mix2, u);

    elseif contains(ProblemType, '_R') && length(range) > 1
        mix3 = mix2;
        mix2 = self.PS.str2;
        self.Misc.config.labely = 'Temperature $T$ [K]';

        if contains(ProblemType, {'DRIVEN'})

            if contains(ProblemType, 'OVER')
                labelx = 'Over-driven';
            elseif contains(ProblemType, 'UNDER')
                labelx = 'Under-driven';
            end

            self.Misc.config.labelx = sprintf('%s ratio $u_1/u_{cj}$', labelx);
            u = cell2vector(mix1, 'drive_factor');
        else
            self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
            u = cell2vector(mix1, 'u');
        end

        ax = plot_figure('u', u, 'T', mix2, 'config', self.Misc.config, 'labelx', self.Misc.config.labelx);
        self.Misc.config.title = 'Reflected';
        self.Misc.config.legend_name = {'Incident', 'Reflected'};
        plot_figure('u', u, 'T', self.PS.strP, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
        % Plot Molar fractions incident state with incident velocity
        self.Misc.config.title = 'Incident';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        ax = displaysweepresults(self, mix2, u);
        set_title(ax, self.Misc.config);
        % Plot Molar fractions reflected state with incident velocity
        self.Misc.config.title = 'Reflected';
        ax = displaysweepresults(self, mix3, u);
        set_title(ax, self.Misc.config);
        % Plot Hugoniot curves
        if strcmp(ProblemType, {'SHOCK_R'})
            ax = plot_hugoniot(self, mix1, mix2);
            ax = plot_hugoniot(self, mix2, mix3, ax);
            set_legends(ax, self.Misc.config.legend_name, 'config', self.Misc.config)
        end
    
    elseif contains(ProblemType, {'OBLIQUE'}) && length(range) > 1
        % Plot molar fractions
        ax = displaysweepresults(self, mix2, range);
        
        % set title
        self.Misc.config.title = 'Strong shock';
        set_title(ax, self.Misc.config);
        % Plot range variable againts temperature
        ax = plot_figure(range_name, range, 'T', mix2, 'config', self.Misc.config);
        if strcmpi(range_name, 'theta')
            % Plot range variable againts temperature
            mix2_2 = self.PS.str2;
            self.Misc.config.legend_name = {'strong shock', 'weak shock'};
            self.Misc.config.legend_location = {'best'};
            plot_figure(range_name, range, 'T', mix2_2, 'config', self.Misc.config, 'ax', ax, 'color', 'auto');
            % Plot molar fractions
            ax = displaysweepresults(self, mix2_2, range);
            % set title
            self.Misc.config.title = 'Weak shock';
            set_title(ax, self.Misc.config);
        end
        
    elseif contains(ProblemType, {'DRIVEN'}) && length(range) > 1
        % Set label x axis
        if contains(ProblemType, 'OVER')
            labelx = 'Over-driven';
        elseif contains(ProblemType, 'UNDER')
            labelx = 'Under-driven';
        end

        self.Misc.config.labelx = sprintf('%s ratio $u_1/u_{cj}$', labelx);
        % Plot drive factor against temperature
        plot_figure('drive_factor', mix1, 'T', mix2, 'config', self.Misc.config, 'labelx', self.Misc.config.labelx);
        % Plot Molar fractions
        plot_molar_fractions(self, mix1, 'drive_factor', 'Xi', 'ydata', mix2);
        % Plot Hugoniot curves
        plot_hugoniot(self, mix1, mix2);
    end

    if strcmp(ProblemType, {'DET'}) && length(range) > 1
        plot_figure('u', mix1, 'T', mix2, 'config', self.Misc.config);
    end

    if strcmp(ProblemType, {'ROCKET'}) && length(range) > 1

        if FLAG_PLOT_RANGE
            ax = plot_figure('phi', mix1, 'cstar', mix2, 'phi', 'config', self.Misc.config);
            %     legend_name = sprintf('$p = %.2f$ [bar]', self.PS.strR.p);
            %     for i = 2:length(self)
            %         legend_name{i} = sprintf('$p = %.2f$ [bar]', self{i}.PS.strR.p);
            %         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'cstar', self.Misc.config, ax, legend_name);
            %     end

            %     self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
            %     self.Misc.config.labely = 'Temperature $T$ [K]';
            %     ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'T', self.Misc.config);
            %     for i = 2:length(self)
            %         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'T', self.Misc.config, ax, legend_name);
            %     end

            self.Misc.config.legend_name = {'$I_{\rm sp}$', '$I_{\rm vac}$'};
            ax = plot_figure('phi', mix1, 'I_sp', mix2, 'config', self.Misc.config);
            plot_figure('phi', mix1, 'I_vac', mix2, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');

            ax = plot_figure('OF', mix1, 'I_sp', mix2, 'config', self.Misc.config);
            plot_figure('OF', mix1, 'I_vac', mix2, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
        end

        if isfield(self.Misc.FLAGS_PROP, 'Aratio')

            if self.Misc.FLAGS_PROP.Aratio
                self.Misc.config.legend_name = {'$I_{\rm sp}$', '$I_{\rm vac}$'};
                self.Misc.config.labelx = 'Area ratio $A_e/A_t$';
                self.Misc.config.labely = 'Specific impulse $I$ [s]';
                ax = plot_figure(mix2, mix2, 'Aratio', 'I_sp', self.Misc.config);
                plot_figure(mix2, mix2, 'Aratio', 'I_vac', self.Misc.config, ax, 'color', 'auto');
            end

        end

        if isfield(self.Misc.FLAGS_PROP, 'pR')

            if self.Misc.FLAGS_PROP.pR
                self.Misc.config.labelx = 'Pressure $p$ [bar]';
                self.Misc.config.labely = 'Characteristic velocity $c^*$ [m/s]';
                plot_figure('p', mix1, 'cstar', mix2, 'config', self.Misc.config);

                self.Misc.config.legend_name = {'$I_{\rm sp}$', '$I_{\rm vac}$'};
                self.Misc.config.labelx = 'Pressure $p$ [bar]';
                self.Misc.config.labely = 'Specific impulse $I$ [s]';
                ax = plot_figure('p', mix1, 'I_sp', mix2, 'config', self.Misc.config);
                plot_figure('p', mix1, 'I_vac', mix2, 'config', self.Misc.config, 'axes', ax, 'color', 'auto');
            end

        end

        %     for i = 2:length(self)
        %         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'I_sp', self.Misc.config, ax, legend_name);
        %     end
    end

    % EXPORT RESULTS
    export_results(self);
end

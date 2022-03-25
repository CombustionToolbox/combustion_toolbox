function postResults(self)
% Abbreviations ---------------------
mix1 = self.PS.strR;
mix2 = self.PS.strP;
phi = self.PD.phi.value;
timer_0 = self.Misc.timer_0;
ProblemType = self.PD.ProblemType;
% -----------------------------------
self.Misc.timer_loop = toc(timer_0);

fprintf('Elapsed time is %.6g seconds\n', self.Misc.timer_loop);

titlename = strrep(ProblemType, '_', ' ');
self.Misc.config.tit = titlename;
FLAG_PLOT_PHI = numel(phi)>1 && all(phi(2:end) ~= phi(1));

if ~isfield(self.Misc.FLAGS_PROP, 'TP')
    self.Misc.FLAGS_PROP.TP = false;
end

% PLOTS REACTANTS
if self.Misc.FLAGS_PROP.TR && length(phi) > 1
    self.Misc.config.labelx = 'Temperature reactants $T$ [K]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(self, mix2, self.PD.range);
elseif self.Misc.FLAGS_PROP.pR && length(phi) > 1
    self.Misc.config.labelx = 'Pressure reactants $p$ [bar]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(self, mix2, self.PD.range);
elseif isfield(self.Misc.FLAGS_PROP, 'pP')
    if self.Misc.FLAGS_PROP.pP && length(phi) > 1
        self.Misc.config.labelx = 'Pressure $p$ [bar]';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        displaysweepresults(self, mix2, self.PD.range);
    end
end

% PLOTS PRODUCTS
if FLAG_PLOT_PHI
    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(self, mix2, phi);
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        self.Misc.config.labely = 'Temperature $T$ [K]';
        plot_figure(phi, mix2,'phi','T',self.Misc.config,self.PD.CompleteOrIncomplete);
    end
end

if self.Misc.FLAGS_PROP.TP && length(phi) > 1
    self.Misc.config.labelx = 'Temperature $T$ [K]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(self, mix2, self.PD.range);

elseif (strcmp(ProblemType,{'HP'}) || strcmp(ProblemType,{'EV'})) && length(phi) > 1
    self.Misc.config.labelx = 'Adiabatic flame temperature $T$ [K]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    T = cell2vector(mix2, 'T');
    displaysweepresults(self, mix2, T);

elseif strcmp(ProblemType,{'SP'}) && length(phi) > 1
    self.Misc.config.labelx = 'Pressure $p$ [bar]';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(self.PD.pP.value, mix2,'p','T',self.Misc.config,self.PD.CompleteOrIncomplete);

elseif strcmp(ProblemType,{'SV'}) && length(phi) > 1
    self.Misc.config.labelx = 'Volume ratio $v_P/v_R$';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(self.PD.vP_vR.value, mix2,'v_P/v_R','T',self.Misc.config,self.PD.CompleteOrIncomplete);  
    T = cell2vector(mix2, 'T');
    displaysweepresults(self, mix2, T);
elseif strcmp(ProblemType,{'SHOCK_POLAR'})
    % Shock polars
    plot_shock_polar(self, mix1, mix2);
%     % Molar fractions [-] againts incident velocity [m/s] 
%     self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
%     self.Misc.config.labely = 'Molar fraction $X_i$';
%     for i = length(mix2):-1:1
%         beta = cell2vector(mix2{i}.polar, 'beta');
%         displaysweepresults(self, mix2{i}.polar, beta);
%     end
elseif strcmp(ProblemType,{'SHOCK_POLAR_R'})
    mix3 = mix2;
    mix2 = self.PS.str2;
    mix2_1 = self.PS.str2_1;
    mix3_1 = self.PS.str3_1;
    mix3_2 = self.PS.str3_2;
    % Shock polars incident
    plot_shock_polar(self, mix1, mix2);
    % Shock polars reflected
    plot_shock_polar(self, mix2, mix3, mix2_1, mix1);
    % molar fractions [-] againts wave angle [deg]
    self.Misc.config.labelx = 'Wave angle $\beta$ [deg]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    for i = length(mix2):-1:1
        titlename = 'Deflection angle $\theta';
        self.Misc.config.tit = sprintf('%s = %.2f$ [deg]', titlename, mix2_1{i}.theta);
        beta = cell2vector(mix3{i}.polar, 'beta');
        ax = displaysweepresults(self, mix3{i}, beta);
        set_title(ax, self.Misc.config);
%         xlim(ax, [min(mix2{i}.polar.beta), 90]);
    end
    % Temperature [K] againts deflection angle [deg]
    titlename = '\mathcal{M}_1';
    legend_name = {'state 2', 'state 3 - weak', 'state 3 - strong'};
    self.Misc.config.tit = sprintf('%s = %.2f', titlename, mix1{i}.u / mix1{i}.sound);
    self.Misc.config.labelx = 'Deflection $\theta$ [deg]';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    ax = plot_figure(mix2_1, mix2_1, 'theta', 'T',self.Misc.config, self.PD.CompleteOrIncomplete); 
    plot_figure(mix3_1, mix3_1, 'theta', 'T',self.Misc.config, self.PD.CompleteOrIncomplete, ax); 
    plot_figure(mix3_2, mix3_2, 'theta', 'T',self.Misc.config, self.PD.CompleteOrIncomplete, ax); 
    set_legends(ax, legend_name, self.Misc.config)
    % Pressure [bar] againts deflection angle [deg]
    self.Misc.config.labely = 'Pressure $p$ [K]';
    ax = plot_figure(mix2_1, mix2_1, 'theta', 'p',self.Misc.config, self.PD.CompleteOrIncomplete); 
    plot_figure(mix3_1, mix3_1, 'theta', 'p',self.Misc.config, self.PD.CompleteOrIncomplete, ax); 
    plot_figure(mix3_2, mix3_2, 'theta', 'p',self.Misc.config, self.PD.CompleteOrIncomplete, ax); 
    set_legends(ax, legend_name, self.Misc.config)
elseif strcmp(ProblemType,{'SHOCK_I'}) && length(phi) > 1
    self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(self.PD.u1.value, mix2,'u','T',self.Misc.config,self.PD.CompleteOrIncomplete);
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

elseif contains(ProblemType, '_R') && length(phi) > 1
    mix3 = mix2;
    mix2 = self.PS.str2;
    self.Misc.config.labely = 'Temperature $T$ [K]';
    if contains(ProblemType, {'OVERDRIVEN'})
        self.Misc.config.labelx = 'Overdriven ratio $u_1/u_{cj}$';
        u = cell2vector(mix1, 'overdriven');
    else
        self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
        u = cell2vector(mix1, 'u');
    end
    ax = plot_figure(u, mix2,'u','T',self.Misc.config,self.PD.CompleteOrIncomplete); 
    self.Misc.config.tit = 'Reflected';
    legend_name = {'Incident', 'Reflected'};
    plot_figure(u,self.PS.strP,'u','T',self.Misc.config,self.PD.CompleteOrIncomplete, ax, legend_name); 
    % Plot Molar fractions incident state with incident velocity
    self.Misc.config.tit = 'Incident';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    ax = displaysweepresults(self, mix2, u);
    set_title(ax, self.Misc.config);
    % Plot Molar fractions reflected state with incident velocity
    self.Misc.config.tit = 'Reflected';
    ax = displaysweepresults(self, mix3, u);
    set_title(ax, self.Misc.config);
    % Plot Hugoniot curves
    if strcmp(ProblemType,{'SHOCK_R'})
        ax = plot_hugoniot(self, mix1, mix2);
        ax = plot_hugoniot(self, mix2, mix3, ax);
        set_legends(ax, legend_name, self.Misc.config)
    end

elseif strcmp(ProblemType,{'DET_OVERDRIVEN'}) && length(phi) > 1
    self.Misc.config.labelx = 'Overdriven ratio $u_1/u_{cj}$';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(mix1, mix2, 'overdriven', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
    % Plot Molar fractions
    self.Misc.config.labelx = 'Overdriven ratio $u_1/u_{cj}$ [m/s]';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    overdriven = cell2vector(mix1, 'overdriven');
    displaysweepresults(self, mix2, overdriven);
    % Plot Hugoniot curves
    plot_hugoniot(self, mix1, mix2);
end

if strcmp(ProblemType,{'DET'}) && length(phi) > 1
    self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(mix1, mix2, 'u', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
end

if strcmp(ProblemType,{'ROCKET'}) && length(phi) > 1
    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Characteristic velocity $c^*$ [m/s]';
    ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'cstar', self.Misc.config, self.PD.CompleteOrIncomplete);
%     legend_name = sprintf('$p = %.2f$ [bar]', self.PS.strR.p); 
%     for i = 2:length(self)
%         legend_name{i} = sprintf('$p = %.2f$ [bar]', self{i}.PS.strR.p); 
%         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'cstar', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
%     end

%     self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
%     self.Misc.config.labely = 'Temperature $T$ [K]';
%     ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
%     for i = 2:length(self)
%         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
%     end

    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Specific impulse $I_{sp}$ [m/s]';
    ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'I_sp', self.Misc.config, self.PD.CompleteOrIncomplete);
%     for i = 2:length(self)
%         ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'I_sp', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
%     end
end
% EXPORT RESULTS
export_results(self);
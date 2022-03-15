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

self.Misc.config.tit = ProblemType;
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
if ~strcmp(ProblemType,'DET_OVERDRIVEN') && FLAG_PLOT_PHI
    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(self, mix2, phi);
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        self.Misc.config.labely = 'Temperature $T$ [K]';
        plot_figure(phi, mix2,'phi','T',self.Misc.config,self.PD.CompleteOrIncomplete);
    end

elseif self.Misc.FLAGS_PROP.TP && length(phi) > 1
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

elseif strcmp(ProblemType,{'SHOCK_R'}) && length(phi) > 1
    mix3 = mix2;
    mix2 = self.PS.str2;
    self.Misc.config.labelx = 'Incident velocity $u_1$ [m/s]';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    self.Misc.config.tit = 'SHOCK_I';
    ax = plot_figure(self.PD.u1.value, mix2,'u','T',self.Misc.config,self.PD.CompleteOrIncomplete); 
    self.Misc.config.tit = 'SHOCK_R';
    plot_figure(self.PD.u1.value,self.PS.strP,'u','T',self.Misc.config,self.PD.CompleteOrIncomplete, ax, {'$SHOCK_I$', '$SHOCK_R$'}); 
    ax = plot_hugoniot(self, mix1, mix2);
    ax = plot_hugoniot(self, mix1, mix3, ax);
    ax = plot_hugoniot(self, mix2, mix3, ax);

elseif strcmp(ProblemType,{'DET_OVERDRIVEN'}) && length(phi) > 1
    self.Misc.config.labelx = 'Overdriven ratio $u_1/u_{cj}$';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    plot_figure(mix1, mix2, 'overdriven', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
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
    legend_name = sprintf('$p = %.2f$ [bar]', self.PS.strR.p); 
    for i = 2:length(self)
        legend_name{i} = sprintf('$p = %.2f$ [bar]', self{i}.PS.strR.p); 
        ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'cstar', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
    end

    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Temperature $T$ [K]';
    ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
    for i = 2:length(self)
        ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
    end

    self.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    self.Misc.config.labely = 'Specific impulse $I_sp$ [m/s]';
    ax = plot_figure(self.PS.strR, self.PS.strP, 'phi', 'I_sp', self.Misc.config, self.PD.CompleteOrIncomplete);
    for i = 2:length(self)
        ax = plot_figure(self{i}.PS.strR, self{i}.PS.strP, 'phi', 'I_sp', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
    end
end
% EXPORT RESULTS
export_results(self);
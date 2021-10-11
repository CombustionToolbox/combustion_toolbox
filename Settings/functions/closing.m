function closing(app)
% Abbreviations ---------------------
mix1 = app.PS.strR;
mix2 = app.PS.strP;
phi = app.PD.phi.value;
display_species = app.Misc.display_species;
timer_0 = app.Misc.timer_0;
LS = app.S.LS;
mintol_display = app.C.mintol_display;
ProblemType = app.PD.ProblemType;
% -----------------------------------

disp('TIME:')
toc(timer_0);

app.Misc.config.tit = ProblemType;

if numel(phi)>1 && all(phi(2:end) ~= phi(1)) && ~strcmp(ProblemType,'DET_OVERDRIVEN')
    app.Misc.config.labelx = 'Equivalence Ratio $\phi$';
    app.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(app, mix2, phi);
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(phi, mix2,'phi','T',app.Misc.config,app.PD.CompleteOrIncomplete);
    end

elseif strcmp(ProblemType,{'TP'}) && length(phi) > 1
    app.Misc.config.labelx = 'Temperature $T [K]$';
    app.Misc.config.labely = 'Molar fraction $X_i$';
    displaysweepresults(app, mix2, app.PD.range);

elseif strcmp(ProblemType,{'SP'}) && length(phi) > 1
    app.Misc.config.labelx = 'Pressure $p [bar]$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(app.PD.pP.value, mix2,'p','T',app.Misc.config,app.PD.CompleteOrIncomplete);

elseif strcmp(ProblemType,{'SV'}) && length(phi) > 1
    app.Misc.config.labelx = 'Volume ratio $v_P/v_R$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(app.PD.vP_vR.value, mix2,'v_P/v_R','T',app.Misc.config,app.PD.CompleteOrIncomplete);  

elseif strcmp(ProblemType,{'SHOCK_I'}) && length(phi) > 1
    app.Misc.config.labelx = 'Incident velocity $u_1$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(app.PD.u1.value, mix2,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete);
    plot_hugoniot(app, mix1, mix2);
        
elseif strcmp(ProblemType,{'SHOCK_R'}) && length(phi) > 1
    mix3 = mix2;
    mix2 = app.PS.str2;
    app.Misc.config.labelx = 'Incident velocity $u_1$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    app.Misc.config.tit = 'SHOCK_I';
    ax = plot_figure(app.PD.u1.value, mix2,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete); 
    app.Misc.config.tit = 'SHOCK_R';
    plot_figure(app.PD.u1.value,app.PS.strP,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete, ax, {'$SHOCK_I$', '$SHOCK_R$'}); 
    ax = plot_hugoniot(app, mix1, mix2);
    ax = plot_hugoniot(app, mix1, mix3, ax);
    ax = plot_hugoniot(app, mix2, mix3, ax);

elseif strcmp(ProblemType,{'DET_OVERDRIVEN'}) && length(phi) > 1
    app.Misc.config.labelx = 'Overdriven ratio $u_1/u_{cj}$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(mix1, mix2, 'overdriven', 'T', app.Misc.config, app.PD.CompleteOrIncomplete);
    plot_hugoniot(app, mix1, mix2);
% elseif numel(phi)>1 && all(phi(2:end) == phi(1))
%     app.Misc.config.labelx = 'Critical Equivalence Ratio $\phi_c$';
%     app.Misc.config.labely = 'Temperature $T [K]$';
%     plot_figure(mix2, mix2,'phi_c','T',app.Misc.config,app.PD.CompleteOrIncomplete);
end

if strcmp(ProblemType,{'DET'}) && length(phi) > 1
    app.Misc.config.labelx = 'Incident velocity $u_1$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(mix1, mix2, 'u', 'T', app.Misc.config, app.PD.CompleteOrIncomplete);
end
%% EXCEL I/O
ExportExcel(app);
function closing(app)
% Abbreviations ---------------------
strP = app.PS.strP;
phi = app.PD.phi.value;
display_species = app.M.display_species;
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
    displaysweepresults(app, strP, phi);
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(phi,strP,'phi','T',app.Misc.config,app.PD.CompleteOrIncomplete);
    end
elseif any(strcmp(ProblemType,{'SP'}))
        app.Misc.config.labelx = 'Pressure $p [bar]$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.pP.value,strP,'p','T',app.Misc.config,app.PD.CompleteOrIncomplete);
elseif any(strcmp(ProblemType,{'SV'}))
        app.Misc.config.labelx = 'Volume ratio $v_P/v_R$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.vP_vR.value,strP,'v_P/v_R','T',app.Misc.config,app.PD.CompleteOrIncomplete);     
elseif any(strcmp(ProblemType,{'SHOCK_I'}))
        app.Misc.config.labelx = 'Incident velocity $u_1$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.u1.value,strP,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete);
        plot_hugoniot(app);
elseif any(strcmp(ProblemType,{'SHOCK_R'}))
        app.Misc.config.labelx = 'Incident velocity $u_1$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        app.Misc.config.tit = 'SHOCK_I';
        ax = plot_figure(app.PD.u1.value,app.PS.str2,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete); 
        app.Misc.config.tit = 'SHOCK_R';
        plot_figure(app.PD.u1.value,app.PS.strP,'u','T',app.Misc.config,app.PD.CompleteOrIncomplete, ax, {'$SHOCK_I$', '$SHOCK_R$'}); 
elseif numel(phi)>1 && all(phi(2:end) == phi(1))
    app.Misc.config.labelx = 'Critical Equivalence Ratio $\phi_c$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(strP,strP,'phi_c','T',app.Misc.config,app.PD.CompleteOrIncomplete);
end
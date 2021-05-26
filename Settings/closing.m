function closing(app,strP,phi,display_species,timer_0,LS,mintol_display,ProblemType)
disp('TIME:')
toc(timer_0);

app.Misc.config.tit = ProblemType;

if numel(phi)>1 && all(phi(2:end) ~= phi(1)) && ~strcmp(ProblemType,'DET_OVERDRIVEN')
    if isempty(display_species) 
        displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$');
    else
        displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$',display_species);
    end
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        app.Misc.config.labelx = 'Equivalence Ratio $\phi [-]$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.phi.value,app.PS.strP,'phi','T',app.Misc.config,app.PD.CompleteOrIncomplete);
    end
elseif any(strcmp(ProblemType,{'SP'}))
        app.Misc.config.labelx = 'Pressure $p [bar]$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.pP.value,app.PS.strP,'p','T',app.Misc.config,app.PD.CompleteOrIncomplete);
elseif any(strcmp(ProblemType,{'SV'}))
        app.Misc.config.labelx = 'Volume ratio $v_P/v_R$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.vP_vR.value,app.PS.strP,'v_P/v_R','T',app.Misc.config,app.PD.CompleteOrIncomplete);        
elseif numel(phi)>1 && all(phi(2:end) == phi(1))
    app.Misc.config.labelx = 'Critical Equivalence Ratio $\phi_c$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(app.PS.strP,app.PS.strP,'phi_c','T',app.Misc.config,app.PD.CompleteOrIncomplete);
end
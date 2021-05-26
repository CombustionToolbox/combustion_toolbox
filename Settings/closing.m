function closing(app,strP,phi,display_species,timer_0,LS,mintol_display,ProblemType)
disp('TIME:')
toc(timer_0);
if numel(phi)>1 && all(phi(2:end) ~= phi(1)) && ~strcmp(ProblemType,'DET_OVERDRIVEN')
    if isempty(display_species) 
        displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$');
    else
        displaysweepresults(strP,phi,LS,mintol_display,'Equivalence Ratio $\phi$',display_species);
    end
    if ~any(strcmp(ProblemType,{'TP','TV'}))
        app.Misc.config.tit = ProblemType;
        app.Misc.config.labelx = 'Equivalence Ratio $\phi [-]$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.phi.value,app.PS.strP,'phi','T',app.Misc.config,app.PD.CompleteOrIncomplete);
    end
elseif numel(phi)>1 && all(phi(2:end) == phi(1))
    app.Misc.config.tit = ProblemType;
    app.Misc.config.labelx = 'Critical Equivalence Ratio $\phi_c$';
    app.Misc.config.labely = 'Temperature $T [K]$';
    plot_figure(app.PS.strP,app.PS.strP,'phi_c','T',app.Misc.config,app.PD.CompleteOrIncomplete);
end
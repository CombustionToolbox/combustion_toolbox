function Plotting(app,display_species,flag,flag_TP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   app = all data
% OUTPUT:
%   Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(app.PD.phi.value)>1 && all(app.PD.phi.value(2:end) ~= app.PD.phi.value(1)) && flag
    if isempty(display_species)
        displaysweepresults(app.PS.strP,app.PD.phi.value,app.S.NameSpecies,app.C.mintol_display,'Equivalence Ratio $\phi$');
    else
        displaysweepresults(app.PS.strP,app.PD.phi.value,app.S.NameSpecies,app.C.mintol_display,'Equivalence Ratio $\phi$',display_species);
    end
    if ~any(strcmp(app.ProblemType.Value,{'1','4'}))
        app.Misc.config.tit = app.Reactants.Items{sscanf(app.Reactants.Value,'%d')};
        app.Misc.config.labelx = 'Equivalence Ratio $\phi [-]$';
        app.Misc.config.labely = 'Temperature $T [K]$';
        plot_figure(app.PD.phi.value,app.PS.strP,'phi','T',app.Misc.config,app.Reaction.Items{sscanf(app.Reaction.Value,'%d')});
%         display_Tp_phi(app.PS.strP,app.PD.phi.value,app.Reactants.Items{sscanf(app.Reactants.Value,'%d')},app.Reaction.Items{sscanf(app.Reaction.Value,'%d')});
    end
end
if flag_TP
    if isempty(display_species)
        displaysweepresults(app.PS.strP,cell2vector(app.PS.strP,'T'),app.S.NameSpecies,app.C.mintol_display,'Temperature $T_P$');
    else
        displaysweepresults(app.PS.strP,cell2vector(app.PS.strP,'T'),app.S.NameSpecies,app.C.mintol_display,'Temperature $T_P$',display_species);
    end
end
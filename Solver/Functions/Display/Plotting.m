function Plotting(app,display_species)
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
if numel(app.PD.phi.Value)>1 && all(app.PD.phi.Value(2:end) ~= phi(1)) && flag
    if isempty(display_species)
        displaysweepresults(app.PS.strP,app.PD.phi.Value,app.S.NameSpecies,app.C.mintol_display);
    else
        displaysweepresults(app.PS.strP,app.PD.phi.Value,app.S.NameSpecies,app.C.mintol_display,display_species);
    end
    if ~any(strcmp(app.ProblemType.Value,{'1','4'}))
        display_Tp_phi(app.PS.strP,app.PD.phi.Value,app.Reactants.Items{sscanf(numberReactants,'%d')},app.Reaction.Items{aux3});
    end
end
end
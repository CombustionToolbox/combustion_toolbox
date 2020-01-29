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

function vector = struct2vector(str,field)
    Nstruct = length(str);
    for i=Nstruct:-1:1
        vector(i) = str{i}.(field);
    end
end
function plot_figure(x,y,config)
config.linewidth;
config.fontsize;
config.tit 
config.leg = {varargin{4}};

%%% CHECK TIT FOR LATEX
config.tit = strrep(config.tit,'%','\%');
%%%
% Plot configuration
set(figure,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize,'BoxStyle','full')
grid on; box on; hold on; axis tight

xlabel(config.labelx,'FontSize',config.fontsize,'interpreter','latex');
ylabel(config.labely,'FontSize',config.fontsize,'interpreter','latex');

plot(x,y,'LineWidth',config.linewidth);
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');
title({config.tit},'Interpreter','latex','FontSize',config.fontsize+4);

% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');
end
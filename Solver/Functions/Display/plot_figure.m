function plot_figure(x,y,config)
config.linewidth;
config.fontsize;
config.tit 
% config.leg = {varargin{4}};

%%% CHECK TIT FOR LATEX
config.tit = strrep(config.tit,'%','\%');
%%%
% Plot configuration
f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize,'BoxStyle','full')
movegui(f,'center');
grid on; box on; hold on; axis tight

xlabel(config.labelx,'FontSize',config.fontsize,'interpreter','latex');
ylabel(config.labely,'FontSize',config.fontsize,'interpreter','latex');

plot(x,y,'LineWidth',config.linewidth);
legend(leg,'FontSize',fontsize,'Location','northeast','interpreter','latex');
title({config.tit},'Interpreter','latex','FontSize',config.fontsize+4);

% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');
end
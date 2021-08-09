function ax = plot_hugoniot(varargin)
app = varargin(1);
app = app{1,1};
app.Misc.config.labelx = '$R^{-1}$';
app.Misc.config.labely = '$P$';
config = app.Misc.config;

        
if nargin > 1
    ax = varargin(2);
    ax = ax{1,1};
else
    ax = set_figure(config);
end
%%% CHECK TIT FOR LATEX
config.tit = strrep(config.tit,'%','\%');
%%%
x1 = cell2vector(app.PS.strR, 'rho');
x2 = cell2vector(app.PS.strP, 'rho');
y1 = cell2vector(app.PS.strR, 'p');
y2 = cell2vector(app.PS.strP, 'p');
x = x2./x1;
y = y2./y1;
% Plot configuration
    plot(ax, 1./x, y, 'LineWidth', config.linewidth);
end

function ax = set_figure(config)
    f = figure;
    set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
        'outerposition',[0.1 0.1 0.9 0.8])
    ax = axes(f);
    set(ax,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')
    grid(ax, 'on'); box(ax, 'on'); hold(ax, 'on'); axis(ax, 'tight');
    xlabel(ax, config.labelx,'FontSize',config.fontsize,'interpreter','latex');
    ylabel(ax, config.labely,'FontSize',config.fontsize,'interpreter','latex');
    title({strcat('$',config.tit,'$')},'Interpreter','latex','FontSize',config.fontsize+4);
    set(ax,'yscale','log')
    xlim(ax, [0, 1])
    ylim(ax, [1, 10000])
end


function ax = plot_hugoniot(varargin)
app  = varargin(1); app = app{1,1};
mix1 = varargin(2); mix1 = mix1{1,1};
mix2 = varargin(3); mix2 = mix2{1,1};

app.Misc.config.labelx = '$R^{-1}$';
app.Misc.config.labely = '$P$';
config = app.Misc.config;

        
if nargin > 3
    ax = varargin(4);
    ax = ax{1,1};
    legend
else
    ax = set_figure(config);
end
%%% CHECK TITLE COMPATIBILITY LATEX
config.tit = strrep(config.tit,'%','\%');
%%%
x1 = cell2vector(mix1, 'rho');
x2 = cell2vector(mix2, 'rho');
y1 = cell2vector(mix1, 'p');
y2 = cell2vector(mix2, 'p');
x = x2./x1;
y = y2./y1;
% Plot configuration
plot(ax, 1./x, y, 'LineWidth', config.linewidth);
end

%%% SUB-PASS FUNCTIONS
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


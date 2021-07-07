function plot_figure(varargin)
if nargin == 4
    if isstruct(varargin{1})
        x = struct2vector(varargin{1},varargin{3});
    else
        x = varargin{1};
    end
    if isstruct(varargin{2})
        y = struct2vector(varargin{2},varargin{4});
    else
        y = varargin{2};
    end
    config.linewidth = 1.8;
    config.fontsize  = 18;
    config.tit = {'Figure'};
    config.labelx = '$x$';
    config.labely = '$y$';
end
if nargin == 5
    if isstruct(varargin{1})
        x = struct2vector(varargin{1},varargin{3});
    else
        x = varargin{1};
    end
    if isstruct(varargin{2})
        y = struct2vector(varargin{2},varargin{4});
    else
        y = varargin{2};
    end
    config = varargin{5};
end
if nargin == 6
    if iscell(varargin{1})
        x = cell2vector(varargin{1},varargin{3});
    else
        x = varargin{1};
    end
    if iscell(varargin{2})
        y = cell2vector(varargin{2},varargin{4});
    else
        y = varargin{2};
    end
    config = varargin{5};
    leg = varargin(6);
end

%%% CHECK TIT FOR LATEX
config.tit = strrep(config.tit,'%','\%');
%%%
% Plot configuration
f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')

grid on; box on; hold on; axis tight

xlabel(config.labelx,'FontSize',config.fontsize,'interpreter','latex');
ylabel(config.labely,'FontSize',config.fontsize,'interpreter','latex');
xlim([min(x),1.02*max(x)])
ylim([min(y),1.02*max(y)])
plot(x,y,'LineWidth',config.linewidth);
if nargin == 4
    legend(leg,'FontSize',config.fontsize-2,'Location','northeast','interpreter','latex');
end
title({strcat('$',config.tit,'$')},'Interpreter','latex','FontSize',config.fontsize+4);

% movegui(f,'center');
% filename2 = strcat(fpath,filename);
% saveas(fig,filename2,'epsc');
end
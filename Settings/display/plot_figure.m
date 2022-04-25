function ax = plot_figure(varargin)

    % Default
    config.linewidth = 1.8;
    config.fontsize  = 18;
    config.title = 'Figure';
    config.labelx = '$x$';
    config.labely = '$y$';

    if nargin == 4
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
    end
    if nargin == 5
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
    end
    if nargin > 5
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
    if nargin > 6
        ax = varargin(7);
        ax = ax{1,1};
    else
        ax = set_figure(config, x, y);
    end
    if nargin == 8
        legend_name = varargin(8);
        legend_name = legend_name{1,1};
    end
    %%% CHECK TIT FOR LATEX
    config.title = strrep(config.title,'%','\%');
    %%%
    % Plot configuration
    plot(ax, x, y, 'LineWidth', config.linewidth);
    if nargin == 8
        legend(ax, legend_name,'FontSize',config.fontsize-2,'Location','northeastoutside','interpreter','latex');
    end
end

% SUB-PASS FUNCTIONS
function ax = set_figure(config, x, y)
    f = figure;
    set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
        'outerposition',[0.1 0.1 0.9 0.8])
    ax = axes(f);
    set(ax,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')
%     grid(ax, 'on'); box(ax, 'on');
    hold(ax, 'on'); axis(ax, 'tight');
    xlabel(ax, config.labelx,'FontSize',config.fontsize,'interpreter','latex');
    ylabel(ax, config.labely,'FontSize',config.fontsize,'interpreter','latex');
    title({strcat('$',config.title,'$')},'Interpreter','latex','FontSize',config.fontsize+4);
% 	xlim(ax, [min(x),1.02*max(x)])
%     ylim(ax, [min(y),1.02*max(y)])
end
function displaysweepresults(self, str, xvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SWEEP RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   self = App class with all the data
%   str = Prop. of state x (phi,species,...)
%   xvar = x variable
% OUTPUT:
%   Plot: Molar fractions of str for the range xvar given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help displaysweepresults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
nfrec = 3;

% Get inputs
display_species = self.Misc.display_species;
NameSpecies = self.S.LS;
mintol = self.C.mintol_display;

% Function
if length(xvar)>1
    axes = set_figure(self);
    [yvar, all_ind, indlabel] = get_parameters(str, NameSpecies, display_species, mintol);
    plot_line(self, axes, xvar, yvar, all_ind, indlabel, NameSpecies)
    set_limits_axes(axes, xvar, yvar, mintol);
    set_legends(self, axes, NameSpecies(all_ind));
end

end
function axes = set_figure(self)
    % Set figure
%     if ~self.Misc.FLAG_GUI
        f = figure;
        set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
            'outerposition',[0.05 0.05 0.9 0.9]);
        axes = gca;
        set(axes,'LineWidth',self.Misc.config.linewidth,'FontSize',self.Misc.config.fontsize-2,'BoxStyle','full')
%     else
%         axes = self.UIAxes;
%         cla(axes);
%     end
    grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
    
    xlabel(axes, self.Misc.config.labelx,'FontSize',self.Misc.config.fontsize,'interpreter','latex');
    ylabel(axes, self.Misc.config.labely,'FontSize',self.Misc.config.fontsize,'interpreter','latex');
end

function plot_line(self, axes, xvar, yvar, indy, indlabel, label_name)
    NE = numel(indy);
    maxLdisplay = self.Misc.config.colorpaletteLenght;
    if NE > maxLdisplay
        NUM_COLORS = maxLdisplay;
    else
        NUM_COLORS = NE;
    end
    LINE_STYLES = {'-', '--', ':', '-.'};
    NUM_STYLES  = length(LINE_STYLES);

    colorbw = brewermap(NUM_COLORS, self.Misc.config.colorpalette);

    k = 1;
    z = 1;
    for i=1:numel(indy)
        dl = plot(axes, xvar, yvar(indy(i),:), 'LineWidth', self.Misc.config.linewidth, 'color', colorbw(k,:), 'LineStyle', LINE_STYLES{z});
        label_line(self, axes, dl, label_name{indlabel(i)}, i)
        k = k + 1;
        if k == maxLdisplay
            k = 1;
            z = z + 1;
            if z > NUM_STYLES
                z = 1;
            end
        end
    end
end

function label_line(self, axes, dl, label_name, i)
    if self.Misc.FLAG_LABELS
        if mod(i,nfrec)==0
            loc_label = 'right';
        else
            loc_label = 'left';
        end
        label(dl, label_name, 'FontSize', self.Misc.config.fontsize,...
            'location', loc_label, 'axes', axes);
    end
end

function set_limits_axes(axes, xvar, yvar, mintol)
    set(axes,'yscale','log')
    xlim(axes, [min(xvar),max(xvar)])
    ymin = 10^floor(log(abs(min(yvar(yvar>0))))/log(10));
    if ymin > mintol
        ylim(axes, [ymin,1])
    else
        ylim(axes, [mintol,1])
    end
end

function set_legends(self, axes, legends_name)
    legend(axes, legends_name,'FontSize', self.Misc.config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
end

function [yvar, all_ind, indlabel] = get_parameters(str, NameSpecies, display_species, mintol)
    all_ind = [];
    Nstruct = length(str);
    yvar = zeros(length(str{1}.Xi), Nstruct);
    if isempty(display_species)
        for i=Nstruct:-1:1
            yvar(:,i) = str{i}.Xi;
            j = str{i}.Xi > mintol;
            ind = find(j>0);
            all_ind = [all_ind; ind];
        end
        all_ind = unique(all_ind);
        indlabel = all_ind;
    else
        for i=Nstruct:-1:1
            yvar(:,i) = str{i}.Xi;
        end
        for i = 1:length(NameSpecies)
            if any(strcmpi(NameSpecies{i}, display_species))
                all_ind = [all_ind, i];
            end
        end
        indlabel = any(yvar(all_ind,:)'>mintol);
        all_ind = all_ind(indlabel);
    end
end

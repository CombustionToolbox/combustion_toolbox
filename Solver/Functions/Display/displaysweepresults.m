function displaysweepresults(app, str, xvar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY SWEEP RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varangin = [str,phi,display_species]
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
display_species = app.M.display_species;
NameSpecies = app.S.LS;
mintol = app.C.mintol_display;

% Function
if length(xvar)>1
    axes = set_figure(app);
    [yvar, all_ind, indlabel] = get_parameters(str, NameSpecies, display_species, mintol);
    plot_line(app, axes, xvar, yvar, all_ind, indlabel, NameSpecies)
    set_limits_axes(app, axes, xvar, yvar, mintol);
    set_legends(app, axes, NameSpecies(all_ind));
end

end
function axes = set_figure(app)
    % Set figure
%     if ~app.Misc.FLAG_GUI
        f = figure;
        set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
            'outerposition',[0.05 0.05 0.9 0.9]);
        axes = gca;
        set(axes,'LineWidth',app.Misc.config.linewidth,'FontSize',app.Misc.config.fontsize-2,'BoxStyle','full')
%     else
%         axes = app.UIAxes;
%         cla(axes);
%     end
    grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
    
    xlabel(axes, app.Misc.config.labelx,'FontSize',app.Misc.config.fontsize,'interpreter','latex');
    ylabel(axes, app.Misc.config.labely,'FontSize',app.Misc.config.fontsize,'interpreter','latex');
end

function plot_line(app, axes, xvar, yvar, indy, indlabel, label_name)
    NE = numel(indy);  
    colorbw = brewermap(NE, 'Spectral');
    for i=1:numel(indy)
        dl = plot(axes, xvar, yvar(indy(i),:), 'LineWidth', app.Misc.config.linewidth, 'color', colorbw(i,:));
        label_line(app, axes, dl, label_name{indlabel(i)}, i)
    end
end

function label_line(app, axes, dl, label_name, i)
    if app.Misc.FLAG_LABELS
        if mod(i,nfrec)==0
            loc_label = 'right';
        else
            loc_label = 'left';
        end
        label(dl, label_name, 'FontSize', app.Misc.config.fontsize,...
            'location', loc_label, 'axes', axes);
    end
end

function set_limits_axes(app, axes, xvar, yvar, mintol)
    set(axes,'yscale','log')
    xlim(axes, [min(xvar),max(xvar)])
    ymin = 10^floor(log(abs(min(yvar(yvar>0))))/log(10));
    if ymin > mintol
        ylim(axes, [ymin,1])
    else
        ylim(axes, [mintol,1])
    end
end

function set_legends(app, axes, legends_name)
    legend(axes, legends_name,'FontSize',app.Misc.config.fontsize-6,'Location','northeastoutside','interpreter','latex');
end

function [yvar, all_ind, indlabel] = get_parameters(str, NameSpecies, display_species, mintol)
    all_ind = [];
    Nstruct = length(str);
    yvar = zeros(length(str.Xi),Nstruct);
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

function plot_molar_fractions_validation(results1, results2, varname_x, varname_y, species, varargin)
    % Default values
    nfrec = 1;
    mintol_display = 1e-14;
    config.fontsize = 18;
    config.linewidth = 1.8;
    
    dataname_x = get_dataname(varname_x);
    dataname_y = get_dataname(varname_y);
    results1.(varname_x) = cell2vector(select_data(results1, dataname_x), varname_x);
    results1.(varname_y) = cell2vector(select_data(results1, dataname_y), varname_y);
    index_species_CT = find_ind(results1.Misc.LS_original, species);

    f = figure;
    set(f,'units','normalized','innerposition',[0.05 0.05 0.9 0.9],...
        'outerposition',[0.05 0.05 0.9 0.9]);
    axes = gca;
    set(axes,'LineWidth',config.linewidth,'FontSize',config.fontsize-2,'BoxStyle','full')
    grid(axes, 'off'); box(axes, 'off'); hold(axes, 'on'); axes.Layer = 'Top';
    xlabel(axes, 'Equivalence ratio, $\phi$','FontSize',config.fontsize,'interpreter','latex');
    ylabel(axes, 'Molar fraction, $X_i$','FontSize',config.fontsize,'interpreter','latex');
    set(axes,'yscale','log')
    xlim(axes, [min(results1.(varname_x)), max(results1.(varname_x))])
    ylim(axes, [mintol_display, 1])
    colorbw = brewermap(length(species), 'Spectral');
    for i=1:length(index_species_CT)
        plot(results1.(varname_x), results1.(varname_y)(index_species_CT(i), :), 'LineWidth', config.linewidth, 'color', colorbw(i,:));
    end
    for i=1:length(species)
        plot(results2.(varname_x)(1:nfrec:end), results2.(varname_y)(i, 1:nfrec:end), 'd', 'LineWidth', config.linewidth, 'color', colorbw(i,:));
    end
    legendname = species;
    legend(legendname, 'FontSize', config.fontsize-6, 'Location', 'northeastoutside', 'interpreter', 'latex');
end

function dataname = get_dataname(var)
    if strcmpi(var, 'phi')
        dataname = 'PD.phi.value';
    else
        dataname = 'PS.strP';
    end
end

function dataselected = select_data(self, dataname)
    index = strfind(dataname, '.');
    index = [index, length(dataname) + 1];
    N = length(index);
    dataselected = self;
    pos1 = 1;
    for i = 1:N
        pos2 = index(i) - 1;
        varname = dataname(pos1:pos2);
        dataselected = dataselected.(varname);
        pos1 = index(i) + 1;
    end
end
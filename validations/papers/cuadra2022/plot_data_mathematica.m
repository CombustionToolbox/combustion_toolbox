filename1 = 'C:\Users\user\Google Drive\Phd\Shocks\Dissociation Turbulence\NumericAir\range_vorticity';
prefix1 = 'dataW3D';
sufix1 = 'zetaair.txt';

filename2 = 'C:\Users\user\Google Drive\Phd\Shocks\Dissociation Turbulence\NumericAir\range_PDF';
prefix2 = 'dataPDF';
sufix2 = 'air.txt';

FLAG_SAVE = true;
% Miscellaneous
greensfmc = [56, 122, 102]/255;
redsfmc = [156, 61, 43]/255;

% Load data
for i = 200:-1:1
    data_1(:, :, i) = importdata(fullfile(filename1, strcat(prefix1, sprintf('%d', i), sufix1)));
    data_2(:, :, i) = importdata(fullfile(filename2, strcat(prefix2, sprintf('%d', i), sufix2)));
end

% Plot results
[ax_1, config, fig] = set_figure;
config.labelx = '$\mathcal{M}_1$';
config.labely = '$|\Omega|^2$';
[ax_1, config] = set_figure(ax_1, config);

config.labelx = '$\mathcal{M}_1$';
config.labely = 'PDF';
ax_2 = axes('Position', [.22 .7 .3 .2]);
ax_2 = set_figure(ax_2, config);

start = 100; % zeta = 1
offset = 0.01;

j_20k = 5105;


set_default_values_1(ax_1);
set_default_values_2(ax_2);

for i = 1:start-1
    if i == 0
        color_long = [0, 0, 0];
        color_short = [0, 0, 0];
    else
        color_long = greensfmc;
        color_short = redsfmc;
    end

    % Clear axes
    cla(ax_1);
    cla(ax_2);

    % Background
    patch(ax_1, [1 2 2 1], [100 100 0 0], [78, 78, 88]/255, 'FaceAlpha', 0.1, 'LineStyle', 'none')
    patch(ax_1, [2 13 13 2], [100 100 0 0], [117, 178, 178]/255, 'FaceAlpha', 0.3, 'LineStyle', 'none')
    patch(ax_1, [13 30 30 13], [100 100 0 0], [198, 174, 104]/255, 'FaceAlpha', 0.3, 'LineStyle', 'none')
    patch(ax_1, [30 57 57 30], [100 100 0 0], [101, 17, 0]/255, 'FaceAlpha', 0.3, 'LineStyle', 'none')

    % Main axes
    plot(ax_1, data_1(1:j_20k-1, 1, start - i), data_1(1:j_20k-1, 2, start - i), 'color', color_long, 'LineWidth', 2.5);
    plot(ax_1, data_1(1:j_20k-1, 1, start + i), data_1(1:j_20k-1, 2, start + i), 'color', color_short, 'LineWidth', 2.5);

    plot(ax_1, data_1(j_20k+1:end, 1, start - i), data_1(j_20k+1:end, 2, start - i), '--', 'color', color_long, 'LineWidth', 2.5);
    plot(ax_1, data_1(j_20k+1:end, 1, start + i), data_1(j_20k+1:end, 2, start + i), '--', 'color', color_short, 'LineWidth', 2.5);
%     plot(ax, data(:, 1, start), data(:, 2, start), 'k', 'LineWidth', 1.2);
    long_label = strcat('\zeta = ', sprintf('%.2f', 1 - offset * i));
    short_label = strcat('\zeta = ', sprintf('%.2f', 1 + offset * i));
    config.title = ['\color[rgb]{0.22,0.49,0.40}', long_label, ',     \color[rgb]{0.61,0.24,0.17}', short_label];
    title(ax_1, config.title, 'Interpreter', 'tex', 'fontsize', config.fontsize + 2);
    
    % Inset
    plot(ax_2, data_2(:, 1, start - i), data_2(:, 2, start - i), 'color', color_long, 'LineWidth', 2.5);
    plot(ax_2, data_2(:, 1, start + i), data_2(:, 2, start + i), 'color', color_short, 'LineWidth', 2.5);
%     plot(ax_2, data_2(1:j_20k-1, 1, start - i), data_2(1:j_20k-1, 2, start - i), 'color', color_long, 'LineWidth', 2.5);
%     plot(ax_2, data_2(1:j_20k-1, 1, start + i), data_2(1:j_20k-1, 2, start + i), 'color', color_short, 'LineWidth', 2.5);
%     plot(ax_2, data_2(j_20k+1:end, 1, start - i), data_2(j_20k+1:end, 2, start - i), '--', 'color', color_long, 'LineWidth', 2.5);
%     plot(ax_2, data_2(j_20k+1:end, 1, start + i), data_2(j_20k+1:end, 2, start + i), '--', 'color', color_short, 'LineWidth', 2.5);


    if FLAG_SAVE
        folderpath = strcat(pwd,'\Validations\papers\cuadra2022\Figures\');
        stack_trace = dbstack;
        filename1 = stack_trace.name;
        saveas(fig, strcat(folderpath, prefix1, sprintf('%d', i)), 'svg');
    end
%     pause(0.1);
end

function set_default_values_1(ax)
    box(ax, 'off');
    xlim(ax, [1, 57]);
    ylim(ax, [0, 100]);
    set(ax, 'xscale', 'log');
%     set(ax, 'yscale', 'log');
    xticks(ax, [1, 5, 10, 20, 40])
    yticks(ax, [0, 20, 40, 60, 80, 100])

    hYLabel = get(ax, 'YLabel');
    set(hYLabel, 'rotation', 0, 'VerticalAlignment', 'middle')
end

function set_default_values_2(ax)
    box(ax, 'off');
    xlim(ax, [1, 57]);
    ylim(ax, [0, 1]);
    set(ax, 'xscale', 'log');
    xticks(ax, [1, 5, 10, 20, 40])
%     yticks(ax, [0, 20, 40, 60, 80, 100])
    hYLabel = get(ax, 'YLabel');
    set(hYLabel, 'rotation', 90, 'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center')
end
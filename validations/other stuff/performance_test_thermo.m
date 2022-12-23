function performance_test_thermo()
    % Performance test thermodynamic polynomials:
    %   * Combustion Toolbox custom polynomials.
    %   * NASA's 9 coefficient polynomials.

    % Inputs
    N = 650;
    properties = {'cp', 'h0', 'g0', 's0'};
    % Initialization
    self = App('NASA ALL');
    LS = self.S.LS;
    NS = self.S.NS;
    NP = length(properties);
    % Compute performance Test
    for i = NP:-1:1
        [time_NASA, time_CT] = compute_performance(self, LS, properties{i}, N, NS);
        sum_time_NASA(i) = sum(time_NASA);
        sum_time_CT(i) = sum(time_CT);
%         print_results(time_NASA, time_CT, N, NS);
%         plot_time_cases(time_NASA, time_CT, NS);
    end
    % Bar plot with the total execution time [s] per property
    plot_bar(properties, sum_time_NASA, sum_time_CT, N, NS, NP);
end

% SUB-PASS FUNCTIONS
function [time_NASA, time_CT] = compute_performance(self, LS, property, N, NS)
    % Compute performance Test

    % Get thermodynamic functions
    [funname_NASA, funname_CT] = set_inputs_thermo_validations(property);
    % Compute performance test
    for j = NS:-1:1
        species = LS{j};
        temperature = self.DB.(species).T;
        temperature = linspace(temperature(1), temperature(end), N);
    
        tic
        for i=1:N
            funname_NASA(species, temperature(i), self.DB);
        end
        time_NASA(j) = toc;
        
        tic
        for i=1:N
            funname_CT(species, temperature(i), self.DB);
        end
        time_CT(j) = toc;
    end
end

function print_results(time_NASA, time_CT, N, NS)
    % Print results
    fprintf('\n\n');
    fprintf('PERFORMANCE TEST\n');
    fprintf('Computations: %g\n', N * NS);
    fprintf('NASA:         %.4f seconds\n', sum(time_NASA));
    fprintf('CT:           %.4f seconds\n', sum(time_CT));
    fprintf('Performance:  %.2f %%\n', (sum(time_NASA) / sum(time_CT) - 1) * 100);
end

function ax = plot_time_cases(time_NASA, time_CT, NS)
    % Plot time per case

    % Compute max value
    threshold = max(max(time_NASA), max(time_CT));
    % Plot results
    cases = 1:1:NS;
    [ax, config] = set_figure();
    plot(ax, cases, time_CT/threshold, 'k', 'LineWidth', config.linewidth);
    plot(ax, cases, time_NASA/threshold, 'k--', 'LineWidth', config.linewidth);
    xlabel(ax, 'Cases')
    ylabel(ax, 'Time*');
    legendname = {'Combustion Toolbox', 'NASA'};
    legend(ax, legendname, 'FontSize', config.fontsize-6, 'Location', ...
        'northeastoutside', 'interpreter', 'latex');
    ylim(ax, [0, 1]);
end

function ax = plot_bar(properties, sum_time_NASA, sum_time_CT, N, NS, NP)
    % Bar plot with the total execution time [s] per property
    [ax, config] = set_figure();
    config.tit = strcat(sprintf('%d', N * NS), ' computations per case');
    width = .4;
    cases = 1:1:NP;
    bar(ax, cases, sum_time_NASA, width,'FaceColor', [0.2 0.2 0.5]);
    bar(ax, cases, sum_time_CT, width,'FaceColor', [0 0.7 0.7]);
    ax.XTick = cases; 
    ax.XTickLabels = properties;
    xlabel(ax, '');
    ylabel(ax, 'Time [s]');
    legendname = {'NASA', 'Combustion Toolbox'};
    legend(ax, legendname, 'FontSize', config.fontsize-6, 'Location', ...
        'northeast', 'interpreter', 'latex');
    title(ax, config.tit, 'Interpreter', 'latex', 'FontSize', config.fontsize+4);

    for i = NP:-1:1 
        label_NASA = sprintf('%.2f', sum_time_NASA(i));
        label_CT = sprintf('%.2f', sum_time_CT(i));
        text(ax, cases(i), sum_time_NASA(i), label_NASA, 'Color', 'w', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', ...
            'top', 'FontSize', config.fontsize-6)
        text(ax, cases(i), sum_time_CT(i), label_CT, 'Color', 'w', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', ...
            'top', 'FontSize', config.fontsize-6)
    end
end
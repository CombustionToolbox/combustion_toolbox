%% CASE HIGH TEMPERATURE
clear;

n = 100;
p = [logspace(0, 3, n), 1050, 1100];
n = n + 2;

T = 273;

LS = {'H2', 'N2', 'O2', 'C2H4'};
Xi = 1;

self = App;
ax = set_figure;

c = [160, 83, 161;...
     170, 201, 75;...
     232, 74, 79;...
     83, 109, 174;...
     243, 234, 81] / 255;

for j = 1:length(LS)
    for i = n:-1:1
        v_molar(i) = eos_PR(self, LS(j), Xi, T, p(i));
    end
    Z = convert_bar_to_Pa(p) .* v_molar ./ (self.C.R0 * T);
    plot(ax, convert_bar_to_atm(p), Z, 'LineWidth', self.Misc.config.linewidth, 'Color', c(j , :));
    legend_label(j) = {species2latex(LS{j})};
end
plot(ax, [0, max(convert_bar_to_atm(p))], [1, 1], 'k--', 'LineWidth', self.Misc.config.linewidth);

xlim(ax, [0, 1000]);
ylim(ax, [0, 2.5]);
xlabel(ax, 'Pressure [atm]')
ylabel(ax, 'Compressibility factor')
legend(ax, legend_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize, 'Location', 'best');
%% CASE HIGH TEMPERATURE
clear;

n = 100;
p = [logspace(0, 3, n), 1050, 1100];
n = n + 2;

T = [173, 298, 873, 10000];

LS = {'NO2'};
Xi = 1;

self = App;
ax = set_figure;

c = [160, 83, 161;...
     170, 201, 75;...
     232, 74, 79;...
     83, 109, 174;...
     243, 234, 81] / 255;

for j = 1:length(T)
    for i = n:-1:1
        v_molar(i) = eos_PR(self, LS, Xi, T(j), p(i));
    end
    Z = convert_bar_to_Pa(p) .* v_molar ./ (self.C.R0 * T(j));
    plot(ax, convert_bar_to_atm(p), Z, 'LineWidth', self.Misc.config.linewidth, 'Color', c(j , :));
end
plot(ax, [0, max(convert_bar_to_atm(p))], [1, 1], 'k--', 'LineWidth', self.Misc.config.linewidth);

legend_label = {'$T = 173$ K', '$T = 298$ K', '$T = 873$ K', '$T = 10,000$ K', 'Ideal gas'};


xlim(ax, [0, 1000]);
ylim(ax, [0.5, 2]);
xlabel(ax, 'Pressure [atm]')
ylabel(ax, 'Compressibility factor')
legend(ax, legend_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize, 'Location', 'best');
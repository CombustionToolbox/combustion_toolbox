%% CASE DIFFERENT SPECIES
clear;

n = 100;
p = [logspace(0, 3, n), 1050, 1100];
n = n + 2;

T = 273;

LS = {'H2', 'N2', 'O2', 'C2H4', 'CO2'};
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
        V(i) = eos_VanderWaals(self, LS(j), Xi, T, p(i));
    end
    Z = convert_bar_to_Pa(p) .* V ./ (self.C.R0 * T);
    plot(ax, convert_bar_to_atm(p), Z, 'LineWidth', self.Misc.config.linewidth, 'Color', c(j , :));
    legend_label(j) = {species2latex(LS{j})};
end
plot(ax, [0, max(convert_bar_to_atm(p))], [1, 1], 'k--', 'LineWidth', self.Misc.config.linewidth);

title_label = sprintf('Deviation of Z with pressure for different species at $T = %g$ K using PR EoS', T);

xlim(ax, [0, 1000]);
ylim(ax, [0, 2.5]);
xlabel(ax, 'Pressure, $p$ [atm]')
ylabel(ax, 'Compressibility factor, $Z$')
legend(ax, legend_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize, 'Location', 'best');
title(ax, title_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize)
%% CASE HIGH TEMPERATURE
clear;

n = 100;
p = [logspace(0, 3, n), 1050, 1100];
n = n + 2;

T = [173, 298, 873, 10000];

LS = {'N2'};
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
        V(i) = eos_VanderWaals(self, LS, Xi, T(j), p(i));
    end
    Z = convert_bar_to_Pa(p) .* V ./ (self.C.R0 * T(j));
    plot(ax, convert_bar_to_atm(p), Z, 'LineWidth', self.Misc.config.linewidth, 'Color', c(j , :));
end
plot(ax, [0, max(convert_bar_to_atm(p))], [1, 1], 'k--', 'LineWidth', self.Misc.config.linewidth);

legend_label = {'$T = 173$ K', '$T = 298$ K', '$T = 873$ K', '$T = 10,000$ K', 'Ideal gas'};


xlim(ax, [0, 1000]);
ylim(ax, [0.5, 2]);
xlabel(ax, 'Pressure, $p$ [atm]')
ylabel(ax, 'Compressibility factor, $Z$')
legend(ax, legend_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize, 'Location', 'best');
title(ax, ['Deviation of Z with temperature and pressure for ', species2latex(LS{1}), ' using PR EoS'], 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize)
%% CASE Van der Waals
clear;

n = 100;
p = [logspace(0, 2, n)];

T = [150, 170, 190.6, 210]; 

LS = {'CH4'};
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
        V(i) = eos_VanderWaals(self, LS, Xi, T(j), p(i));
    end
    v = V;
    plot(ax, v, convert_bar_to_Pa(p), 'LineWidth', self.Misc.config.linewidth, 'Color', c(j , :));
end

legend_label = {'$T = 150$ K', '$T = 170$ K', '$T = 190.6$ K', '$T = 210$ K'};
set(ax, 'Xscale', 'log');

ylim(ax, [0, 1e7]);
xlim(ax, [1e-5, 1e-2]);
ylabel(ax, 'Pressure, $p$ [Pa]')
xlabel(ax, 'Volume, $v$ [m$^3$]')
legend(ax, legend_label, 'Interpreter', 'latex', 'fontsize', self.Misc.config.fontsize, 'Location', 'best');
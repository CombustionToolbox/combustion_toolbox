% EXAMPLE: COMPARISON ROCKET FAC MODELS WITH FROZEN (POST-COMBUSTION) AND
% EQUILIBRIUM CHEMISTRY
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Apr 10 2023
% -------------------------------------------------------------------------

% Definitions
range_sub = 1.001:0.05:3;
range_sup = 1.001:0.1:25.1;

% Initialization
self = App();

%% Case 1: frozen chemistry (post-combustion)
FLAG_FROZEN = true;
% Subsonic region
FLAG_SUBSONIC = true;
[T_frozen_sub, M_frozen_sub] = compute_rocket_FAC(self, FLAG_FROZEN, FLAG_SUBSONIC, range_sub, range_sup);
% Supersonic region
FLAG_SUBSONIC = false;
[T_frozen_sup, M_frozen_sup] = compute_rocket_FAC(self, FLAG_FROZEN, FLAG_SUBSONIC, range_sub, range_sup);

%% Case 2: chemical equilibrium
FLAG_FROZEN = false;
% Subsonic region
FLAG_SUBSONIC = true;
[T_equil_sub, M_equil_sub] = compute_rocket_FAC(self, FLAG_FROZEN, FLAG_SUBSONIC, range_sub, range_sup);
% Supersonic region
FLAG_SUBSONIC = false;
[T_equil_sup, M_equil_sup] = compute_rocket_FAC(self, FLAG_FROZEN, FLAG_SUBSONIC, range_sub, range_sup);

%% Plot results
ax = set_figure();
xlim([0, max(range_sup)]);

% Temperature (left axis) 
yyaxis left
set(ax, 'YColor', [0, 0, 0]);

[ax, dl1] = plot_figure('Aratio', range_sub, 'T', T_frozen_sub, 'ax', ax, 'linestyle', '--', 'color', [0, 0, 0]);
ax = plot_figure('Aratio', range_sup, 'T', T_frozen_sup, 'ax', ax, 'linestyle', '--', 'color', [0, 0, 0]);
[ax, dl2] = plot_figure('Aratio', range_sub, 'T', T_equil_sub, 'ax', ax, 'linestyle', '-', 'color', [0, 0, 0]);
ax = plot_figure('Aratio', range_sup, 'T', T_equil_sup, 'ax', ax, 'linestyle', '-', 'color', [0, 0, 0]);

ax = plot_figure('Aratio', [1 1], 'T', ax.YLim, 'ax', ax, 'linestyle', ':', 'color', [0.2, 0.2, 0.2]);

dl3 = plot(ax, max(range_sub), T_equil_sub(end), 'd', 'MarkerSize', 8, 'LineWidth', self.Misc.config.linewidth, 'color', [0, 0, 0], 'MarkerFaceColor', 'auto');
dl4 = plot(ax, 1, 0.5 * (T_equil_sub(1) + (T_equil_sup(1))), 'o', 'MarkerSize', 8, 'LineWidth', self.Misc.config.linewidth, 'color', [0, 0, 0], 'MarkerFaceColor', 'auto');
plot(ax, 1, 0.5 * (T_frozen_sub(1) + (T_frozen_sup(1))), 'o', 'MarkerSize', 8, 'LineWidth', self.Misc.config.linewidth, 'color', [0, 0, 0], 'MarkerFaceColor', 'auto');

% Pre-shock Mach number (right axis)
yyaxis right
set(ax, 'YColor', self.Misc.config.colorline);

ax = plot_figure('Aratio', range_sub, 'M', M_frozen_sub, 'ax', ax, 'linestyle', '--', 'color', self.Misc.config.colorline);
ax = plot_figure('Aratio', range_sup, 'M', M_frozen_sup, 'ax', ax, 'linestyle', '--', 'color', self.Misc.config.colorline);
ax = plot_figure('Aratio', range_sub, 'M', M_equil_sub, 'ax', ax, 'linestyle', '-', 'color', self.Misc.config.colorline);
ax = plot_figure('Aratio', range_sup, 'M', M_equil_sup, 'ax', ax, 'linestyle', '-', 'color', self.Misc.config.colorline);

plot(ax, max(range_sub), M_equil_sub(end), 'd', 'MarkerSize', 8, 'LineWidth', self.Misc.config.linewidth, 'color', self.Misc.config.colorline, 'MarkerFaceColor', 'auto');
plot(ax, 1, 1, 'o', 'MarkerSize', 8, 'LineWidth', self.Misc.config.linewidth, 'color', self.Misc.config.colorline, 'MarkerFaceColor', 'auto');

ax.XLabel.String = 'Area$_i$ / Area throat';

% SUB-PASS FUNCTIONS
function [T, M] = compute_rocket_FAC(self, FLAG_FROZEN, FLAG_SUBSONIC, range_sub, range_sup)
    % Initialize
    self = App('fast', self.DB_master, self.DB, 'Soot formation extended');

    % Initial conditions
    self = set_prop(self, 'TR', 298.15, 'pR', 100 * 1.01325, 'phi', 1.5);
    self.PD.S_Fuel     = {'RP_1'};
    self.PD.S_Oxidizer = {'O2bLb'};
    self.PD.FLAG_IAC = false;
    self.PD.FLAG_FROZEN = FLAG_FROZEN;
    self.PD.FLAG_SUBSONIC = FLAG_SUBSONIC;

    % Additional inputs (depends of the problem selected)
    if FLAG_SUBSONIC
        self = set_prop(self, 'Aratio', range_sub, 'Aratio_c', 3);
    else
        self = set_prop(self, 'Aratio', range_sup, 'Aratio_c', 3);
    end
    
    % Solve problem
    self = solve_problem(self, 'ROCKET');

    % Get results
    mix = self.PS.strP;
    T = cell2vector(mix, 'T');
    u = cell2vector(mix, 'u');
    a = cell2vector(mix, 'sound');
    M = u ./ a;
end
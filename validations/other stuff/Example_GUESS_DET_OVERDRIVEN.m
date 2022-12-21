%% OBTAIN SOLUTIONS AT CHEMICAL EQUILIBRIUM
% Initialize
self = App('Soot Formation');
% Initial conditions
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
% Additional inpunts (depends of the problem selected)
drive_factor = 1:0.1:5;
self = set_prop(self, 'drive_factor', drive_factor);
% Solve problem
self = solve_problem(self, 'DET_OVERDRIVEN');
% Get results
T2 = cell2vector(self.PS.strP, 'T');
p2 = cell2vector(self.PS.strP, 'p');

%% COMPUTE INITIAL GUESSES AS A FUNCTION OF THE DEGREE OF OVERDRIVE
mix2 = [];
for i = length(drive_factor):-1:1
    [p2_guess(i, 3), T2_guess(i, 3)] = guess_shock(self, self.PS.strR{i});
    [p2_guess(i, 2), T2_guess(i, 2)] = guess_det_CEA(self, self.PS.strR{i});
    [p2_guess(i, 1), T2_guess(i, 1)] = guess_det(self, self.PS.strR{i}, drive_factor(i));
end

%% PLOT RESULTS
% Temperature
ax = set_figure; ylim(ax, [2000, 12000])
plot_figure('$\eta$', drive_factor, '$T_2$', T2, 'linestyle', 'k', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$T_2$', T2_guess(:, 1), 'linestyle', '--', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$T_2$', T2_guess(:, 2), 'linestyle', '--', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$T_2$', T2_guess(:, 3), 'linestyle', '--', 'color', 'auto', 'ax', ax);
legend(ax, {'Solution', 'Guess DET CJ - Q', 'Guess DET CJ - CEA', 'Guess shock'}, 'Interpreter', 'latex')
% Pressure
ax = set_figure; ylim(ax, [0, 900]);
plot_figure('$\eta$', drive_factor, '$p_2$', p2, 'linestyle', 'k', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$p_2$', p2_guess(:, 1), 'linestyle', '--', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$p_2$', p2_guess(:, 2), 'linestyle', '--', 'color', 'auto', 'ax', ax);
plot_figure('$\eta$', drive_factor, '$p_2$', p2_guess(:, 3), 'linestyle', '--', 'color', 'auto', 'ax', ax);
legend(ax, {'Solution', 'Guess DET CJ - Q', 'Guess DET CJ - CEA', 'Guess shock'}, 'Interpreter', 'latex')

% NESTED FUNCTIONS
function [p2, T2] = guess_det(self, mix1, drive_factor)
    % Estimate post-shock state considering h2 = h1 + u1^2 / 2
    [p2p1, T2T1] = det_compute_guess(self, mix1, mix1.phi, drive_factor);
    T2 = T2T1 * mix1.T; % [K]
    p2 = p2p1 * mix1.p; % [bar]
end

function [p2, T2] = guess_det_CEA(self, mix1)
    % Estimate post-shock state considering h2 = h1 + u1^2 / 2
    [p2p1, T2T1] = det_compute_guess_CEA(self, mix1);
    T2 = T2T1 * mix1.T; % [K]
    p2 = p2p1 * mix1.p; % [bar]
end

function [p2, T2] = guess_shock(self, mix1)
    % Estimate post-shock state considering h2 = h1 + u1^2 / 2
    M1 = mix1.u / mix1.sound;
    p2p1 = (2 * mix1.gamma * M1^2 - mix1.gamma + 1) / (mix1.gamma + 1);
    mix1.h = (enthalpy_mass(mix1) * 1e3 + velocity_relative(mix1)^2/2) * mass(mix1) * 1e-3; % [kJ]
    self.PD.ProblemType = 'HP';
    mix2 = equilibrate(self, mix1, p2p1 * mix1.p);
    T2 = mix2.T; % [K]
    p2 = p2p1 * mix1.p; % [bar]
end
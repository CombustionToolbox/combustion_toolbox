% -------------------------------------------------------------------------
% EXAMPLE: DET_POLAR_SONIC_AND_MAX
%
% Compute detonation polar curves at standard conditions, a set of 30 species
% considered, and a set of pre-shock Mach numbers M1 = [1.01*M1cj:15], or
% what is the same, drive factors M1/M1_cj = [1.01:2.9069]
%    
% Hydrogen == {'H2O','H2','O2','N2','He','Ar','H','HNO',...
%              'HNO3','NH','NH2OH','NO3','N2H2','N2O3','N3','OH',...
%              'HNO2','N','NH3','NO2','N2O','N2H4','N2O5','O','O3',...
%              'HO2','NH2','H2O2','N3H','NH2NO2'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Oct 07 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Hydrogen');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
self.PD.S_Fuel = {'H2'};
self.PD.S_Oxidizer = {'N2', 'O2'};
self.PD.ratio_oxidizers_O2 = [79, 21]/21;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'drive_factor', 1.01:0.01:2.9069);
%% TUNING PROPERTIES
self.TN.N_points_polar = 300;  % Number of points to compute shock polar
%% SOLVE PROBLEM
self = solve_problem(self, 'DET_POLAR');
%% GET POINTS
theta_sonic = cell2vector(self.PS.strP, 'theta_sonic');
beta_sonic = cell2vector(self.PS.strP, 'beta_sonic');
theta_max = cell2vector(self.PS.strP, 'theta_max');
beta_max = cell2vector(self.PS.strP, 'beta_max');

for i = length(self.PS.strP):-1:1
    postshock_u2n(:, i) = self.PS.strP{i}.polar.un;
    postshock_u2(:, i) = self.PS.strP{i}.polar.u;
    postshock_a(:, i) = self.PS.strP{i}.polar.sound;
    beta(:, i) = self.PS.strP{i}.polar.beta;
    theta(:, i) = self.PS.strP{i}.polar.theta;
end

M2 = postshock_u2 ./ postshock_a;
M2n = postshock_u2n ./ postshock_a;

for i = length(self.PS.strP):-1:1
    ind_row_sonic(i) =  find(M2(:, i) >= 1, 1);
    ind_row_sonic_n(i) =  find(M2n(:, i) >= 1, 1);
    beta_sonic_2(i) = beta(ind_row_sonic(i), i);
    theta_sonic_2(i) = theta(ind_row_sonic(i), i);
end
%% DISPLAY RESULTS (PLOTS)
post_results(self);
%% SMOOTH RESULTS - FOURIER
start_point_sonic = [0 0 0 0 0 0.0508682282375078];
start_point_sonic_2 = [0 0 0 0 0 0.0661385531723257];
start_point_max = [0 0 0 0 0 0.0656724073958934];
[beta_sonic_smooth, theta_sonic_smooth] = smooth_data([90, beta_sonic], [0, theta_sonic], start_point_sonic);
[theta_max_smooth, beta_max_smooth] = smooth_data([0, theta_max], [90, beta_max], start_point_max);
[theta_sonic_2_smooth, beta_sonic_2_smooth] = smooth_data([0, theta_sonic_2], [90, beta_sonic_2], start_point_sonic_2);
%% PLOT
f = figure(2); ax = gca;
plot_figure('$\theta_{\rm M2n sonic}$', theta_sonic_smooth, '$\beta_{\rm M2n sonic}$', beta_sonic_smooth, 'linestyle', 'k:', 'color', 'auto', 'ax', ax)
plot_figure('$\theta_{\rm max}$', theta_max_smooth, '$\beta_{\rm max}$', beta_max_smooth, 'linestyle', 'k:', 'color', 'auto', 'ax', ax)
plot_figure('$\theta_{\rm M2 sonic}$', theta_sonic_2_smooth, '$\beta_{\rm M2 sonic}$', beta_sonic_2_smooth, 'linestyle', 'k:', 'color', 'auto', 'ax', ax)
% -------------------------------------------------------------------------
% EXAMPLE: HP - EFFECT OF PRESSURE
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1 bar) for lean to rich natural gas-air mixtures at
% standard conditions, a set of equivalence ratios phi contained in
% (0.5, 3.5) [-], and two different sets of species "Complete" or "Incomplete"
% The incomplete combustion case is evaluated at different pressures to see
% Le Chatelier's principle, i.e., an increase of pressure shifts
% equilibrium to the side of the reaction with fewer number of moles of
% gas (higher pressures implie less dissociation). A comparison of the
% adiabatic flame temperature, adiabatic index, and Gibbs free energy is
% included.
% 
% * Complete:
%       - lean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'};                (equivalence ratio < 1)
%       - rich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'};          (equivalence ratio > 1)
%       - soot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb', 'CO2', 'H2O'}; (equivalence ratio > equivalence ratio soot)   
%
% * Incomplete:
%       - Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                            'C2','C2H4','CH','CH3','CH4','CN','H',...
%                            'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 29 2022
% -------------------------------------------------------------------------

%% COMPLETE COMBUSTION
% Initialization
self = App('Complete');
% Set fuel composition 
self.PD.S_Fuel = {'CH4', 'C2H6', 'C3H8', 'C4H10_isobutane', 'H2'};
self.PD.N_Fuel = [0.8, 0.05, 0.05, 0.05, 0.05];
% Set oxidizer composition
self = set_air(self, false);
% Set temperature, pressure and equivalence ratio
self = set_prop(self, 'TR', 300, 'pR', 1, 'phi', 0.5:0.005:3.5);
% Set constant pressure for products
self = set_prop(self, 'pP', self.PD.pR.value);
% Solve Problem
self = solve_problem(self, 'HP');
% Save results
results_complete = self;
%% INCOMPLETE COMBUSTION
pressure_vector = self.PD.pR.value * [1, 10, 100];
for i = length(pressure_vector):-1:1
    % Set product species at equilibrium
    self = App('copy', self, 'Soot formation');
    % Set pressure
    self = set_prop(self, 'pP', pressure_vector(i), 'pP', pressure_vector(i));
    % Solve Problem
    self = solve_problem(self, 'HP');
    % Save results
    results_incomplete{i} = self;
end
%% Comparison of adiabatic flame temperature
mix1 = results_complete.PS.strR; % Is the same as the incomplete case (same initial mixture)
mix2_complete = results_complete.PS.strP; 

phi = cell2vector(mix1, 'phi');

self.Misc.config.title = '\rm{Complete\ vs\ Incomplete\ at\ different\ pressures}';
self.Misc.config.labelx = 'Equivalence ratio $\phi$';
self.Misc.config.labely = 'Temperature $T$ [K]';
legend_name = {'Complete'};

ax = plot_figure('phi', phi, 'T', mix2_complete, 'config', self.Misc.config);

for i = 1:length(pressure_vector)
    mix2_incomplete = results_incomplete{i}.PS.strP; 
    ax = plot_figure('phi', phi, 'T', mix2_incomplete, 'config', self.Misc.config, 'ax', ax, 'color', 'auto');
    legend_name(i+1) = {sprintf('$p_{%d} = %.4g$ bar', i, pressure(mix2_incomplete{1}))};
end

set_legends(ax, legend_name, 'FLAG_SPECIES', false)
ax.Legend.Location = 'northeast';
%% Comparison of the adiabatic index
self.Misc.config.labely = 'Adibatic index $\gamma_s$';

ax = plot_figure('phi', phi, 'gamma_s', mix2_complete, 'config', self.Misc.config);

for i = 1:length(pressure_vector)
    mix2_incomplete = results_incomplete{i}.PS.strP; 
    ax = plot_figure('phi', phi, 'gamma_s', mix2_incomplete, 'config', self.Misc.config, 'ax', ax, 'color', 'auto');
end

set_legends(ax, legend_name, 'FLAG_SPECIES', false)
ax.Legend.Location = 'best';
%% Comparison of the Gibbs energy
self.Misc.config.labely = 'Gibbs free energy $g$ [kJ/kg]';

ax = plot_figure('phi', phi, 'g', mix2_complete, 'config', self.Misc.config);

for i = 1:length(pressure_vector)
    mix2_incomplete = results_incomplete{i}.PS.strP; 
    ax = plot_figure('phi', phi, 'g', mix2_incomplete, 'config', self.Misc.config, 'ax', ax, 'color', 'auto');
end

set_legends(ax, legend_name, 'FLAG_SPECIES', false)
ax.Legend.Location = 'best';
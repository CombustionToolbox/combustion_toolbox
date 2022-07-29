% -------------------------------------------------------------------------
% EXAMPLE: HP COMPLETE VS INCOMPLETE
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich CH4-air mixtures at
% standard conditions, a set of equivalence ratios phi contained in
% (0.5, 3) [-], and two different sets of species "Complete" or "Incomplete"
% 
% * Complete:
%       - lean = {'CO2', 'H2O', 'N2', 'Ar', 'O2'};                (equivalence ratio < 1)
%       - rich = {'CO2', 'H2O', 'N2', 'Ar', 'CO', 'H2'};          (equivalence ratio > 1)
%       - soot = {'N2', 'Ar', 'CO', 'H2', 'Cbgrb', 'CO2', 'H2O'}; (equivalence ratio > equivalence ratio soot)   
%
% * Incomplete:
%       - Soot formation == {'CO2','CO','H2O','H2','O2','N2','He','Ar','Cbgrb',...
%                            'C2','C2H4','CH','CH','CH3','CH4','CN','H',...
%                            'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
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
self = set_prop(self, 'TR', 300, 'pR', 1.01325, 'phi', 0.5:0.02:3);
% Set constant pressure for products
self = set_prop(self, 'pP', self.PD.pR.value);
% Solve Problem
self = SolveProblem(self, 'HP');
% Save results
results_complete = self;
%% INCOMPLETE COMBUSTION
% Set product species at equilibrium
self = App('copy', self, 'Soot formation');
% Solve Problem
self = SolveProblem(self, 'HP');
% Save results
results_incomplete = self;
%% COMPARE RESULTS
mix1 = results_complete.PS.strR; % Is the same as the incomplete case (same initial mixture)
mix2_complete = results_complete.PS.strP; 
mix2_incomplete = results_incomplete.PS.strP;

phi = cell2vector(mix1, 'phi');

self.Misc.config.title = '\rm{Complete\ vs\ Incomplete}';
self.Misc.config.labelx = 'Equivalence ratio $\phi$';
self.Misc.config.labely = 'Temperature $T$ [K]';
legend_name = {'Complete', 'Incomplete'};

ax = plot_figure(phi, mix2_complete, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete);
ax = plot_figure(phi, mix2_incomplete, 'phi', 'T', self.Misc.config, self.PD.CompleteOrIncomplete, ax, legend_name);
ax.Legend.Location = 'northeast';
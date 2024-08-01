% -------------------------------------------------------------------------
% EXAMPLE: HP COMPLETE VS INCOMPLETE
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich natural gas-air mixtures at
% standard conditions, a set of equivalence ratios phi contained in
% (0.5, 3) [-], and two different sets of species "Complete" or "Incomplete"
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
% Last update Jul 25 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Initialize solver
solver = EquilibriumSolver('problemType', 'HP');

%% COMPLETE COMBUSTION

% Define chemical system
system = ChemicalSystem(DB, 'complete');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'CH4', 'C2H6', 'C3H8', 'C4H10_isobutane', 'H2'}, 'fuel', [0.8, 0.05, 0.05, 0.05, 0.05]);
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'equivalenceRatio', 0.5:0.02:3);

% Solve problem
solver.solveArray(mixArray1);

%% INCOMPLETE COMBUSTION

% Define chemical system
system = ChemicalSystem(DB, 'soot formation');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'CH4', 'C2H6', 'C3H8', 'C4H10_isobutane', 'H2'}, 'fuel', [0.8, 0.05, 0.05, 0.05, 0.05]);
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray2 = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'equivalenceRatio', 0.5:0.02:3);

% Solve problem
solver.solveArray(mixArray2);

%% COMPARE RESULTS
ax = solver.report(mixArray1, mixArray2);
legend(ax.Children(end), {'Complete', 'Incomplete'}, 'Interpreter', 'latex', 'FontSize', ax.Children(end).FontSize);

% Another possibility is call directly the next functions:
% ax = plotProperties(repmat({mixArray1(1).rangeName}, 1, 9), mixArray1, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray1, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []});
% ax = plotProperties(repmat({mixArray1(1).rangeName}, 1, 9), mixArray1, {'T', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray2, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}, 'ax', ax);
% legend(ax.Children(end), {'Complete', 'Incomplete'}, 'Interpreter', 'latex', 'FontSize', ax.Children(end).FontSize);
% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute rocket propellant performance considering an Infinite-Area-Chamber
% for lean to rich MHF3-MON25 mixtures at 101.325 bar, a set of 24 species
% considered, a set of equivalence ratios phi contained in (2, 5) [-], and
% an exit area ratio A_exit/A_throat = 3 
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Jul 30 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.rocket.*
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'soot Formation');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'CH6N2bLb', 'N2H4bLb'}, 'fuel', [86, 14], 'weightPercentage');
set(mix, {'N2O4bLb', 'N2O3'}, 'oxidizer', [36.67, 63.33], 'weightPercentage');

% Define properties
mixArray1 = setProperties(mix, 'temperature', 298.15, 'pressure', 100 * 1.01325, 'equivalenceRatio', 2:0.01:5, 'areaRatio', 3);

% Initialize solver
solver = RocketSolver('problemType', 'ROCKET_IAC');

% Solve problem
[mixArray1, mixArray2, mixArray3, mixArray4] = solver.solveArray(mixArray1);

% Plot properties
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray4, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray4, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []});
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray3, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray3, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []}, 'ax', ax1);
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 10), mixArray2, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u'}, mixArray2, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], []}, 'ax', ax1);
leg = legend(ax1.Children(end), {'Exit', 'Throat', 'Chamber'}, 'Interpreter', 'latex', 'FontSize', ax1.Children(end).FontSize);

% Plot molar fractions
plotComposition(mixArray2(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2, 'title', 'Chamber');
plotComposition(mixArray3(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray3, 'title', 'Throat');
plotComposition(mixArray4(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray4, 'title', 'Exit');
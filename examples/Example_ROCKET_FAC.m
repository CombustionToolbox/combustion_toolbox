% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Finite-Area-Chamber (FAC)
%
% Compute rocket propellant performance considering an Finite-Area-Chamber
% with an area ratio of the combustion chamber A_chamber/A_throat = 2
% for lean to rich LH2-LOX mixtures at 101.325 bar, a set of 24 species
% considered, a set of equivalence ratios phi contained in (2, 5) [-], and
% an exit area ratio A_exit/A_throat = 3 
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
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
system = ChemicalSystem(DB, 'hydrogen_L');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'H2bLb'}, 'fuel', 1);
set(mix, {'O2bLb'}, 'oxidizer', 1);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 90, 'pressure', 100 * 1.01325, 'equivalenceRatio', 1:0.05:5, 'areaRatioChamber', 2, 'areaRatio', 3);

% Initialize solver
solver = RocketSolver('problemType', 'ROCKET_FAC');

% Solve problem
[mixArray1, mixArray2_inj, mixArray2_c, mixArray3, mixArray4] = solver.solveArray(mixArray1);

% Plot properties
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray4, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray4, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []});
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray3, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray3, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []}, 'ax', ax1);
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 12), mixArray2_c, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound', 'u', 'I_sp', 'I_vac'}, mixArray2_c, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], [], [], []}, 'ax', ax1);
ax1 = plotProperties(repmat({'equivalenceRatio'}, 1, 10), mixArray2_inj, {'T', 'p', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}, mixArray2_inj, 'basis', {[], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], [], []}, 'ax', ax1);
leg = legend(ax1.Children(end), {'Exit', 'Throat', 'Chamber', 'Injector'}, 'Interpreter', 'latex', 'FontSize', ax1.Children(end).FontSize);

% Plot molar fractions
plotComposition(mixArray2_inj(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2_inj, 'title', 'Injector');
plotComposition(mixArray2_c(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2_c, 'title', 'Chamber');
plotComposition(mixArray3(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray3, 'title', 'Throat');
plotComposition(mixArray4(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray4, 'title', 'Exit');
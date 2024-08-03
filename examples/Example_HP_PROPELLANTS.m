% -------------------------------------------------------------------------
% EXAMPLE: HP PROPELLANTS
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (0.5, 5) [-]
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
% Last update April 02 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'HYDROGEN_L');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'H2bLb'}, 'fuel', 1);
set(mix, {'O2bLb'}, 'oxidizer', 1);

% Define properties
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'equivalenceRatio', 0.2:0.05:5);

% Initialize solver
solver = EquilibriumSolver('problemType', 'HP');

% Solve problem
solver.solveArray(mixArray);

% Plot adiabatic flame temperature
plotFigure('phi', [mixArray.equivalenceRatio], 'T', [mixArray.T]);

% Plot molar fractions
plotComposition(mixArray(1), mixArray, 'equivalenceRatio', 'Xi', 'mintol', 1e-14);
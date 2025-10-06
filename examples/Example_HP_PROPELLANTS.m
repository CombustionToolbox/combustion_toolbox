% -------------------------------------------------------------------------
% EXAMPLE: HP PROPELLANTS
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (p = 1.01325 bar) for lean to rich LH2-LOX mixtures at T = 300 K,
% and a set of equivalence ratios phi contained in (0.5, 5) [-]
%
% See wiki or setListspecies method from ChemicalSystem class for predefined
% sets of species
%
% @author: Alberto Cuadra Lara
%                 
% Last update October 06 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

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

% Generate report
report(solver, mixArray);
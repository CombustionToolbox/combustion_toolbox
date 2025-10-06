% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute rocket propellant performance considering an Infinite-Area-Chamber
% for lean to rich LH2-LOX mixtures at 101.325 bar, a set of equivalence
% ratios phi contained in (1, 5) [-], and an exit area ratio A_exit/A_throat = 3
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
import combustiontoolbox.rocket.*

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
mixArray1 = setProperties(mix, 'temperature', 90, 'pressure', 100 * 1.01325, 'equivalenceRatio', 1:0.05:5, 'areaRatio', 3);

% Initialize solver
solver = RocketSolver('problemType', 'ROCKET_IAC');

% Solve problem
[mixArray1, mixArray2, mixArray3, mixArray4] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2, mixArray3, mixArray4);
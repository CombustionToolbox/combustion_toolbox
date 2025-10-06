% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering an Infinite-Area-Chamber (IAC)
%
% Compute rocket propellant performance considering an Infinite-Area-Chamber
% for lean to rich MHF3-MON25 mixtures at 101.325 bar, a set of equivalence
% ratios phi contained in (2, 5) [-], and an exit area ratio A_exit/A_throat = 3
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
set(mix, {'CH6N2bLb', 'N2H4bLb'}, 'fuel', [86, 14], 'weightPercentage');
set(mix, {'N2O4bLb', 'N2O3'}, 'oxidizer', [36.67, 63.33], 'weightPercentage');

% Define properties
mixArray1 = setProperties(mix, 'temperature', 298.15, 'pressure', 100 * 1.01325, 'equivalenceRatio', 2:0.01:5, 'areaRatio', 3);

% Initialize solver
solver = RocketSolver('problemType', 'ROCKET_IAC');

% Solve problem
[mixArray1, mixArray2, mixArray3, mixArray4] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2, mixArray3, mixArray4);
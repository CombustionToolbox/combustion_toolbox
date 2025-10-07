% -------------------------------------------------------------------------
% EXAMPLE: DET_POLAR
%
% Compute detonation polar curves for a stoichiometric H2-air mixture at 
% standard conditions (T = 300 K, p = 1 atm), and a set of pre-shock Mach
% numbers M1 = [5, 6, 7, 8, 10], or what is the same, drive factors 
% M1/M1_cj = [1.0382, 1.2458, 1.4534, 1.6611, 2.0763]
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
import combustiontoolbox.shockdetonation.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'H2'}, 'fuel', 1);
set(mix, {'N2', 'O2'}, 'oxidizer', [79, 21] / 21);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', [1.0382, 1.2458, 1.4534, 1.6611, 2.0763]);

% Initialize solver
solver = DetonationSolver('problemType', 'DET_POLAR', 'numPointsPolar', 300);

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2);
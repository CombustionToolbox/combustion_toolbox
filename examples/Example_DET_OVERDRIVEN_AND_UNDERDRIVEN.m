% -------------------------------------------------------------------------
% EXAMPLE: DET_OVERDRIVEN_AND_UNDERDRIVEN
%
% Compute pre-shock and post-shock state for a planar underdriven to 
% overdriven detonation considering Chapman-Jouguet (CJ) theory for a
% stoichiometric CH4-air mixture at standard conditions (T = 300 K, p = 1 atm),
% and a set of overdrives contained in (1,10) [-].
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
set(mix, {'CH4'}, 'fuel', 1);
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, 'oxidizer', [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1_over = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);
mixArray1_under = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);

% Initialize solver
solver1 = DetonationSolver('problemType', 'DET_OVERDRIVEN');
solver2 = DetonationSolver('problemType', 'DET_UNDERDRIVEN');

% Solve problem
[mixArray1_over, mixArray2_over] = solver1.solveArray(mixArray1_over);
[mixArray1_under, mixArray2_under] = solver2.solveArray(mixArray1_under);

% Generate reports
report(solver1, mixArray1_over, mixArray2_over);
report(solver2, mixArray1_under, mixArray2_under);
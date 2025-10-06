% -------------------------------------------------------------------------
% EXAMPLE: DET_UNDERDRIVEN
%
% Compute pre-shock and post-shock state for a planar under-driven detonation
% considering Chapman-Jouguet (CJ) theory for a stoichiometric CH4-air
% mixture at standard conditions (T = 300 K, p = 1 atm), and a set of 
% underdrives contained in (1,10) [-].
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);

% Initialize solver
solver = DetonationSolver('problemType', 'DET_UNDERDRIVEN');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2);
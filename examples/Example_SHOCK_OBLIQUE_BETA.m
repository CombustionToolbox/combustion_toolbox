% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_OBLIQUE_BETA
%
% Compute pre-shock and post-shock state for a oblique incident shock wave
% at standard conditions (T1 = 300 K, p1 = 1 atm), a pre-shock Mach number
% M1 = 10, and a set of wave angles beta = [20:1:85] [deg]
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
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'M1', 10, 'beta', 20:1:85);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_OBLIQUE');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2);
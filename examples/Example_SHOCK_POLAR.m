% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR
%
% Compute shock polar plots at standard conditions (T = 300 K, p = 1 atm),
% and a set of initial shock front Mach numbers = [2, 3, 5, 14]
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
system.FLAG_ION = true;

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'M1', [2, 3, 5, 14]);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_POLAR');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2);
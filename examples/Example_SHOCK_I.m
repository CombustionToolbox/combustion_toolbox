% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_I
%
% Compute pre-shock and post-shock state for a planar incident shock wave
% at standard conditions (T1 = 300 K, p1 = 1 atm), and a set of pre-shock
% velocities contained in (400, 12000) [m/s]
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'u1', 400:100:12000);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_I');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2);
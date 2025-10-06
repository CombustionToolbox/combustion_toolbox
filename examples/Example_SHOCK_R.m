% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_R
%
% Compute pre-shock and post-shock state for a planar reflected shock wave
% at standard conditions (T = 300 K, p = 1 atm), and a set of initial shock
% front velocities (u1) contained in (360, 9000) [m/s]
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
u1 = logspace(2, 5, 500); u1 = u1(u1<9000); u1 = u1(u1>=360);
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'u1', u1);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_R');

% Solve problem
[mixArray1, mixArray2, mixArray3] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2, mixArray3);
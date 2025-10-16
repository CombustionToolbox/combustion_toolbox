% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_PRANDTL_MEYER
%
% Compute pre-shock and post-shock state for a Prandtl-Meyer expansion wave
% at temperature T1 = 3000 K, pressure p1 = 1 atm, a free-stream Mach
% number M1 = 1, and a set of wave angles theta = [0:100] [deg]
%
% @author: Alberto Cuadra Lara
%                 
% Last update October 15 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Definitions
FLAG_FROZEN = false;

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79, 21]/21);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 3000, 'pressure', 1, 'M1', 1, 'theta', 100);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_PRANDTL_MEYER', 'FLAG_FROZEN', FLAG_FROZEN, 'numPointsPrandtlMeyer', 300);

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Generate report from mixArray
report(solver, mixArray1, mixArray2)

% Generate report from polar of last mixArray2
mixArray2(end).polar(1).rangeName = 'theta';
report(solver, mixArray2(end).polar, mixArray2(end).polar)
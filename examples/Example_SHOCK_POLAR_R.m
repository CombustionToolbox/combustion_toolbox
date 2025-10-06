% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED
%
% Compute shock polar plots at T1 = 226.65 K and p1 = 0.0117 bar
% (altitude 30 km), an initial shock front Mach number = 20, and
% deflection angle theta = 35
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
set(mix, {'N2', 'O2', 'Ar'}, [78, 21, 1] / 21);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 226.65, 'pressure', 0.0117, 'M1', 20, 'theta', 35);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_POLAR_R');

% Solve problem
[mixArray1, mixArray2, mixArray2_1, mixArray3, mixArray3_1, mixArray3_2] = solver.solveArray(mixArray1);

% Generate report
report(solver, mixArray1, mixArray2, mixArray2_1, mixArray3, mixArray3_1, mixArray3_2);
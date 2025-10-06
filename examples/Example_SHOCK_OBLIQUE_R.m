% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_OBLIQUE_R
%
% Compute pre-shock and post-shock state (incident and reflected) for a
% oblique incident shock wave at standard conditions, a set of 51 species
% considered, a pre-shock Mach number = 10, and a deflection angle
% theta = 20 [deg]
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'M1', 10, 'theta', 20);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_OBLIQUE_R');

% Solve problem
[mixArray1, mixArray2, mixArray3_1, mixArray3_2] = solver.solveArray(mixArray1);
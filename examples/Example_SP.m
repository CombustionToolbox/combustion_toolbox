% -------------------------------------------------------------------------
% EXAMPLE: SP
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of pressure (p = 1:100 atm) for a rich CH4-air mixture
% at an initial temperature T = 300 K, and a equivalence ratio phi 1.5 [-]
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
import combustiontoolbox.equilibrium.*

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
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1.01325 * logspace(0, 2, 200), 'equivalenceRatio', 1.5);

% Initialize solver
solver = EquilibriumSolver('problemType', 'SP');

% Solve problem
solver.solveArray(mixArray);

% Generate report
report(solver, mixArray);
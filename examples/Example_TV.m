% -------------------------------------------------------------------------
% EXAMPLE: TV
%
% Compute equilibrium composition at defined temperature (T = 3000 K) and
% defined specific volume (1 m3/kg) for lean to rich CH4-air mixtures defined
% by a set of equivalence ratios (phi) contained in (0.5, 5) [-]
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
mixArray = setProperties(mix, 'temperature', 3000, 'volume', 1, 'equivalenceRatio', 0.5:0.01:5);

% Initialize solver
solver = EquilibriumSolver('problemType', 'TV');

% Solve problem
solver.solveArray(mixArray);

% Generate report
report(solver, mixArray);
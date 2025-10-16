% -------------------------------------------------------------------------
% EXAMPLE: SV
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of specific volume (vSpecific = 0.5:2 m3/kg) for a lean 
% CH4-air mixture at defined specific entropy, and an equivalence ratio
% phi 0.5 [-]
%   
% See wiki or setListspecies method from ChemicalSystem class for predefined
% sets of species
%
% @author: Alberto Cuadra Lara
%                 
% Last update October 16 2025
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
mixArray = setProperties(mix, 'entropySpecific', mix.sSpecific, 'volume', 0.5:0.01:2, 'equivalenceRatio', 0.5);

% Initialize solver
solver = EquilibriumSolver('problemType', 'SV');

% Solve problem
solver.solveArray(mixArray);

% Generate report
report(solver, mixArray);
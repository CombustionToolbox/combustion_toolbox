% -------------------------------------------------------------------------
% EXAMPLE: EXOPLANET WASP-43b - METALLICITY 1
%
% Compute equilibrium vertical composition with a metallicity 1 of WASP-43b
%   
% URL RESULTS TEA:
% https://github.com/dzesmin/RRC-BlecicEtal-2015a-ApJS-TEA/tree/master/Fig6/WASP43b-solar
%
% @author: Alberto Cuadra Lara
%                 
% Last update April 02 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.SolarAbundances
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

% Definitions
listSpecies = {'C2H2_acetylene', 'C2H4', 'C', 'CH4', 'CO2', 'CO', 'H2', 'H2O', 'H2S', 'H', 'HCN', 'He', 'HS_M', 'N2', 'N', 'NH3', 'O', 'S'};
species = {'H', 'He', 'C', 'N', 'O', 'S'};
metallicity = 1;

% Get initial composition from solar abundances
DB_solar = SolarAbundances();
moles = DB_solar.abundances2moles(species, metallicity);

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, listSpecies);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, species, moles);

% Define properties
mixArray = setProperties(mix, 'temperature', linspace(100, 4000, 300), 'pressure', logspace(-5, 2, 300));

% Initialize solver
solver = EquilibriumSolver('problemType', 'TP');

% Solve problem
solver.solveArray(mixArray);

% Plot molar fractions
plotComposition(mixArray(1), mixArray, 'Xi', 'p', 'mintol', 1e-14, 'ydir', 'reverse', 'xscale', 'log');
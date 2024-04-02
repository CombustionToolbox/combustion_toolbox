% -------------------------------------------------------------------------
% EXAMPLE: TV
%
% Compute equilibrium composition at defined temperature (e.g., 3000 K) and
% defined specific volume for lean to rich CH4-air mixtures at standard
% conditions, a set of 24 species considered and a set of equivalence
% ratios (phi) contained in (0.5, 5) [-]
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or list_species() for more predefined sets of species
%  
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update April 02 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'soot formation');

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
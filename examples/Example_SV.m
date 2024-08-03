% -------------------------------------------------------------------------
% EXAMPLE: SV
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of specific volume (0.5, 2) for a lean CH4-air mixture at
% 700 K, a set of 24 species considered, and a equivalence
% ratio phi 0.5 [-]
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more predefined sets of species
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
import combustiontoolbox.utils.display.*

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
mixArray = setProperties(mix, 'temperature', 700, 'volume', 0.5:0.01:2, 'equivalenceRatio', 0.5);

% Initialize solver
solver = EquilibriumSolver('problemType', 'SV');

% Solve problem
solver.solveArray(mixArray);

% Plot adiabatic flame temperature
plotFigure('vSpecific', mixArray, 'T', mixArray);

% Plot molar fractions
plotComposition(mixArray(1), mixArray, 'vSpecific', 'Xi', 'mintol', 1e-14);
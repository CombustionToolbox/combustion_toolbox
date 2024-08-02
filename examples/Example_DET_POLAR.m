% -------------------------------------------------------------------------
% EXAMPLE: DET_POLAR
%
% Compute detonation polar curves at standard conditions, a set of 30 species
% considered, and a set of pre-shock Mach numbers M1 = [5, 6, 7, 8, 10], or
% what is the same, drive factors M1/M1_cj = [1.0382, 1.2458, 1.4534,...
% 1.6611, 2.0763]
%    
% Hydrogen == {'H2O','H2','O2','N2','He','Ar','H','HNO',...
%              'HNO3','NH','NH2OH','NO3','N2H2','N2O3','N3','OH',...
%              'HNO2','N','NH3','NO2','N2O','N2H4','N2O5','O','O3',...
%              'HO2','NH2','H2O2','N3H','NH2NO2'}
%
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%
% Last update Aug 02 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'hydrogen');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'H2'}, 'fuel', 1);
set(mix, {'N2', 'O2'}, 'oxidizer', [79, 21] / 21);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', [1.0382, 1.2458, 1.4534, 1.6611, 2.0763]);

% Initialize solver
solver = DetonationSolver('problemType', 'DET_POLAR', 'numPointsPolar', 300);

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Plot polars
[ax1, ax2, ax3] = plotPolar(mixArray1, mixArray2);
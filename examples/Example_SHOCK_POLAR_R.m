% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED_THERMO
%
% Compute shock polar plots at T1 = 226.65 K and p1 = 0.0117 bar
% (altitude 30 km), a set of 28 species considered, an initial
% shock front Mach number = 20, deflection angle theta = 35, and different
% thermochemical models (chemical equilibrium and frozen chemistry).
%    
% LS== {'eminus', 'Ar', 'Arplus', 'N', ...
%       'Nplus', 'Nminus', 'NO', 'NOplus', 'NO2', ...
%       'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
%       'N2minus', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
%       'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
%       'O2minus', 'O3'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
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
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Definitions
listSpecies = {'eminus', 'Ar', 'Arplus', 'N', ...
      'Nplus', 'Nminus', 'NO', 'NOplus', 'NO2', ...
      'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
      'N2minus', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
      'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
      'O2minus', 'O3'};

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, listSpecies);

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

% Plot polars - incident
plotPolar(mixArray1, mixArray2);

% Plot polars - reflected
plotPolar(mixArray2_1, mixArray3, mixArray2_1, mixArray1);
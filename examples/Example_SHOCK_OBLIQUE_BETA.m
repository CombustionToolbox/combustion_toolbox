% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_OBLIQUE_BETA
%
% Compute pre-shock and post-shock state for a oblique incident shock wave
% at standard conditions, a set of 51 species considered, a pre-shock Mach 
% number M1 = 10, and a set of wave angles beta = [20:1:85] [deg]
%    
% Air_ions == {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', ...
%              'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', ...
%              'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', ...
%              'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'C3O2', 'N', ...
%              'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', ...
%              'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
%              'N2minus', 'NCN', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
%              'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
%              'O2minus', 'O3'}
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

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'M1', 10, 'beta', 20:1:85);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_OBLIQUE');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Plot Hugoniot curve
plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log');

% Plot molar fractions
plotComposition(mixArray2(1), mixArray1, 'beta', 'Xi', 'mintol', 1e-3, 'y_var', mixArray2);

% Plot properties
properties = {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'};
propertiesBasis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};
plotProperties(repmat({'beta'}, 1, length(properties)), mixArray2, properties, mixArray2, 'basis', propertiesBasis);
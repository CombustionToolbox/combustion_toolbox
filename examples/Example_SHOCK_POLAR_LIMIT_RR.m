% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_LIMIT_RR
%
% Compute limit regular reflections considering a calorically imperfect gas
% at an altitude of 30 km from sea level. The calculations are carried out 
% for a set of pre-shock Mach numbers M1 = [1.75:0.1:20]
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
%                 
% Last update Jun 11 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air_ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 226.65, 'pressure', 0.0117, 'M1', 1.75:0.1:20.45);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_POLAR_LIMITRR', 'tol0', 1e-7);

% Solve problem
[mixArray1, mixArray2, mixArray2_1, mixArray3] = solver.solveArray(mixArray1);

% Plot wave angle [deg] (limit regular reflection) against pre-shock Mach number
ax = plotFigure('Pre-shock Mach number', [mixArray1.mach], 'Wave angle limit RR [deg]', [mixArray2_1.beta], 'linestyle', 'k-');

% Plot deflection angle [deg] (limit regular reflection) against pre-shock Mach number
ax = plotFigure('Pre-shock Mach number', [mixArray1.mach], 'Deflection angle limit RR [deg]', [mixArray2_1.theta], 'linestyle', 'k-');

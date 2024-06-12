% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_I
%
% Compute pre-shock and post-shock state for a planar incident shock wave
% at standard conditions, a set of 16 species considered and a set of
% pre-shock velocities contained in (400, 12000)
%    
% Air == {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4',...
%         'N3','C','CO','CO2','Ar'}
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
system = ChemicalSystem(DB, 'air');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Define properties
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'u1', 400:100:12000);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_I');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Plot Hugoniot curve
plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log');

% Plot post-shock temperature
plotFigure('u', mixArray2, 'T', mixArray2);

% Plot molar fractions
plotComposition(mixArray2(1), mixArray1, 'u', 'Xi', 'mintol', 1e-3, 'y_var', mixArray2);
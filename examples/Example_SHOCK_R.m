% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_R
%
% Compute pre-shock and post-shock state for a planar reflected shock wave
% at standard conditions, a set of 16 species considered and a set of
% initial shock front velocities (u1) contained in (360, 9000) [m/s]
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
u1 = logspace(2, 5, 500); u1 = u1(u1<9000); u1 = u1(u1>=360);
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'u1', u1);

% Initialize solver
solver = ShockSolver('problemType', 'SHOCK_R');

% Solve problem
[mixArray1, mixArray2, mixArray3] = solver.solveArray(mixArray1);

% Plot Hugoniot curve
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '-');
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray3.rho], 'p_2 / p_1', [mixArray3.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '--', 'ax', ax1);
legend(ax1, {'Incident', 'Reflected'}, 'Interpreter', 'latex', 'FontSize', ax1.FontSize);

% Plot post-shock temperature
ax2 = plotFigure('u1', [mixArray1.u], 'T', [mixArray2.T], 'color', [0 0 0], 'linestyle', '-');
ax2 = plotFigure('u1', [mixArray1.u], 'T', [mixArray3.T], 'color', [0 0 0], 'linestyle', '--', 'ax', ax2);
legend(ax2, {'Incident', 'Reflected'}, 'Interpreter', 'latex', 'FontSize', ax2.FontSize);

% Plot molar fractions
plotComposition(mixArray3(1), mixArray1, 'u', 'Xi', 'mintol', 1e-3, 'y_var', mixArray3);
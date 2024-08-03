% -------------------------------------------------------------------------
% EXAMPLE: DET UNDERDRIVEN REFLECTED
%
% Compute pre-shock and post-shock state for a reflected planar underdriven
% detonation considering Chapman-Jouguet (CJ) theory for a stoichiometric 
% CH4-air mixture at standard conditions, a set of 24 species considered
% and a set of overdrives contained in (1,10) [-].
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 28 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);

% Initialize solver
solver = DetonationSolver('problemType', 'DET_UNDERDRIVEN_R');

% Solve problem
[mixArray1, mixArray2, mixArray3] = solver.solveArray(mixArray1);

% Plot Hugoniot curves
ax1 = plotFigure('\rho_1 / \rho_i', [mixArray1.rho] ./ [mixArray2.rho], 'p_i / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '-');
ax1 = plotFigure('\rho_1 / \rho_i', [mixArray1.rho] ./ [mixArray3.rho], 'p_i / p_1', [mixArray3.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '--', 'ax', ax1);
legend(ax1, {'Incident', 'Reflected'}, 'Interpreter', 'latex', 'FontSize', ax1.FontSize);

% Plot post-shock temperature
ax2 = plotFigure('driveFactor', mixArray1, 'T', mixArray2, 'linestyle', '-');
ax2 = plotFigure('driveFactor', mixArray1, 'T', mixArray3, 'linestyle', '--', 'ax', ax2);
legend(ax2, {'Incident', 'Reflected'}, 'Interpreter', 'latex', 'FontSize', ax2.FontSize);

% Plot molar fractions
plotComposition(mixArray2(1), mixArray1, 'driveFactor', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2);
plotComposition(mixArray3(1), mixArray1, 'driveFactor', 'Xi', 'mintol', 1e-14, 'y_var', mixArray3);
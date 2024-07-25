% -------------------------------------------------------------------------
% EXAMPLE: DET_OVERDRIVEN_AND_UNDERDRIVEN
%
% Compute pre-shock and post-shock state for a planar underdriven to 
% overdriven detonation considering Chapman-Jouguet (CJ) theory for a
% stoichiometric CH4-air mixture at standard conditions, a set of 24
% species considered and a set of overdrives contained in (1,10) [-].
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
% Last update Oct 12 2022
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
mixArray1_over = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);
mixArray1_under = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 1:0.1:10);

% Initialize solver
solver1 = DetonationSolver('problemType', 'DET_OVERDRIVEN');
solver2 = DetonationSolver('problemType', 'DET_UNDERDRIVEN');

% Solve problem
[mixArray1_over, mixArray2_over] = solver1.solveArray(mixArray1_over);
[mixArray1_under, mixArray2_under] = solver2.solveArray(mixArray1_under);

% Plot Hugoniot curve
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1_over.rho] ./ [mixArray2_over.rho], 'p_2 / p_1', [mixArray2_over.p] ./ [mixArray1_over.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '-');
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1_under.rho] ./ [mixArray2_under.rho], 'p_2 / p_1', [mixArray2_under.p] ./ [mixArray1_under.p], 'xScale', 'log', 'yScale', 'log', 'linestyle', '--', 'ax', ax1);
legend(ax1, {'Overdriven', 'Underdriven'}, 'Interpreter', 'latex', 'FontSize', ax1.FontSize);

% Plot post-shock temperature
ax2 = plotFigure('driveFactor', mixArray1_over, 'T', mixArray2_over, 'linestyle', '-');
ax2 = plotFigure('driveFactor', mixArray1_under, 'T', mixArray2_under, 'linestyle', '--', 'ax', ax2);
legend(ax2, {'Overdriven', 'Underdriven'}, 'Interpreter', 'latex', 'FontSize', ax1.FontSize);

% Plot molar fractions
plotComposition(mixArray2_over(1), mixArray1_over, 'driveFactor', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2_over);
plotComposition(mixArray2_under(1), mixArray1_under, 'driveFactor', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2_under);
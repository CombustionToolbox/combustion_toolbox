% -------------------------------------------------------------------------
% EXAMPLE: DET_OBLIQUE_THETA
%
% Compute pre-shock and post-shock state for a oblique detonation
% considering Chapman-Jouguet (CJ) theory for a stoichiometric CH4-air
% mixture at standard conditions, a set of 24 species considered, an 
% overdrive of 4 and a set of deflection angles [15:5:50] [deg].
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
% Last update Jul 28 2024
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 1, 'driveFactor', 4, 'theta', 15:5:50);

% Initialize solver
solver = DetonationSolver('problemType', 'DET_OBLIQUE');

% Solve problem
[mixArray1, mixArray2_1, mixArray2_2] = solver.solveArray(mixArray1);

% Plot Hugoniot curves
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2_1.rho], 'p_2 / p_1', [mixArray2_1.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'color', 'auto');
ax1 = plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2_2.rho], 'p_2 / p_1', [mixArray2_2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log', 'color', 'auto', 'ax', ax1);
legend(ax1, {'Weak detonation', 'Strong detonation'}, 'Interpreter', 'latex', 'FontSize', ax1.FontSize);

% Plot molar fractions
plotComposition(mixArray2_1(1), mixArray1, 'theta', 'Xi', 'mintol', 1e-3, 'y_var', mixArray2_1);
plotComposition(mixArray2_2(1), mixArray1, 'theta', 'Xi', 'mintol', 1e-3, 'y_var', mixArray2_2);

% Plot properties
properties = {'T', 'p', 'rho', 'h', 'e', 'g', 's', 'gamma_s'};
propertiesBasis = {[], [], [], 'mi', 'mi', 'mi', 'mi', []};
ax2 = plotProperties(repmat({'theta'}, 1, length(properties)), mixArray2_1, properties, mixArray2_1, 'basis', propertiesBasis);
ax2 = plotProperties(repmat({'theta'}, 1, length(properties)), mixArray2_2, properties, mixArray2_2, 'basis', propertiesBasis, 'ax', ax2);
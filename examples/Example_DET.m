% -------------------------------------------------------------------------
% EXAMPLE: DET
%
% Compute pre-shock and post-shock state for a planar detonation
% considering Chapman-Jouguet (CJ) theory for lean to rich CH4-air mixtures
% at standard conditions, a set of 24 species considered and a set of
% equivalence ratios (phi) contained in (0.5, 5) [-]
%   
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
% Last update Jul 25 2024
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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'equivalenceRatio', 0.5:0.01:4);

% Initialize solver
solver = DetonationSolver('problemType', 'DET');

% Solve problem
[mixArray1, mixArray2] = solver.solveArray(mixArray1);

% Plot Hugoniot curve
plotFigure('\rho_1 / \rho_2', [mixArray1.rho] ./ [mixArray2.rho], 'p_2 / p_1', [mixArray2.p] ./ [mixArray1.p], 'xScale', 'log', 'yScale', 'log');

% Plot post-shock temperature
plotFigure('equivalenceRatio', mixArray1, 'T', mixArray2);

% Plot molar fractions
plotComposition(mixArray2(1), mixArray1, 'equivalenceRatio', 'Xi', 'mintol', 1e-14, 'y_var', mixArray2);
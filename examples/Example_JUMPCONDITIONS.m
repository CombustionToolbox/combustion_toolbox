% -------------------------------------------------------------------------
% EXAMPLE: JUMPCONDITIONS
%
% Compute jump conditions for a planar incident shock wave at temperature
% 300 K and pressure 1.01325 bar for an air mixture with 79% N2 and 21% O2
% in % volume and a range of pre-shock Mach numbers M1 = [1, 8]
% 
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Dec 16 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Definitions
FLAG_PAPER = true;
caloricGasModel = CaloricGasModel.imperfect;

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79, 21] / 21);

% Define properties
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1.01325, 'mach', 1:0.1:8);

% Initialize solver
solver = JumpConditionsSolver('caloricGasModel', caloricGasModel, 'FLAG_PAPER', FLAG_PAPER, 'FLAG_RESULTS', false);

% Solve problem
jumpData = solver.solve(mixArray);

% Generate report
report(solver, jumpData);
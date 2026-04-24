% -------------------------------------------------------------------------
% EXAMPLE: JUMPCONDITIONS_PENG_ROBINSON
%
% Compute jump conditions for a planar incident shock wave using both the
% ideal gas and Peng-Robinson equations of state. The initial state is
% biphenyl (C12H10_biphenyl) at T = 746.63996 K and p = 11.637175 bar, with
% pre-shock Mach numbers in the range M1 = [1, 5].
%
% @author: Alberto Cuadra Lara
%          Universidad Carlos III de Madrid
%
% Last update Apr 10 2026
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% User definitions
temperature = 746.63996; % [K]
pressure = 1163717.5e-5; % [bar]
mach = 1:0.01:5;
species = {'C12H10_biphenyl'};
moles = 1;

% Definitions
FLAG_PAPER = true;
caloricGasModel = CaloricGasModel.thermallyPerfect;

% Get Nasa database
DB = NasaDatabase();

% Add Peng-Robinson parameters
DB.species = addPengRobinsonProperties(DB.species);

% Define chemical system
system = ChemicalSystem(DB);

% Define equation of state
eosIdeal = EquationStateIdealGas();
eosPR = EquationStatePengRobinson();

% Initialize mixture
mixIdeal = Mixture(system, 'eos', eosIdeal);
mixPR = Mixture(system, 'eos', eosPR);

% Define chemical state
set(mixIdeal, species, moles);
set(mixPR, species, moles);

% Define properties
mixArrayIdeal = setProperties(mixIdeal, 'temperature', temperature, 'pressure', pressure, 'mach', mach);
mixArrayPR = setProperties(mixPR, 'temperature', temperature, 'pressure', pressure, 'mach', mach);

% Initialize solver
solver = JumpConditionsSolver('caloricGasModel', caloricGasModel, 'FLAG_PAPER', FLAG_PAPER, 'FLAG_RESULTS', false);

% Solve problem
jumpDataIdeal = solver.solve(mixArrayIdeal);
jumpDataPR = solver.solve(mixArrayPR);

% Generate report
report(solver, jumpDataIdeal);
report(solver, jumpDataPR);
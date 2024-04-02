% Example
% Import the necessary classes ( included in edit(fullfile(userpath,'startup.m')) )
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.ShockSolver

tic

% Definitions
M1 = linspace(1.01, 10, 1000);
N = length(M1);

% Load NASA's database
DB = NasaDatabase();

% Define the chemical system, i.e., the set of species that may be present in the chemical transformation at equilibrium
system = ChemicalSystem(DB, 'air_ions');

% Define the initial mixture composition
mix1 = Mixture(system);
set(mix1, {'N2', 'O2'}, [79/21, 1]);

% Define state
mix1Array = mix1.setProperties('temperature', 300, 'pressure', 1, 'mach', 20, 'theta', 35);

% Create the equilibrium solver
solver = ShockSolver('problemType', 'SHOCK_POLAR_R', 'flag_results', true);

% Solve shock incident
[mix1Array, mix2Array, mix3Array, mix4Array] = solver.solveArray(mix1Array);

toc

ax = plot_figure('theta', [mix2Array.theta], 'P', [mix2Array.p] ./ [mix1Array.p]);
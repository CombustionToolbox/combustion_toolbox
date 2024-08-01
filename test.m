% Example
% Import the necessary classes ( included in edit(fullfile(userpath,'startup.m')) )
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

tic

% Definitions
temperature = linspace(1000, 6000, 1000);
N = length(temperature);

% Load NASA's database
DB = NasaDatabase();

% Define the chemical system, i.e., the set of species that may be present in the chemical transformation at equilibrium
system = ChemicalSystem(DB, 'Soot formation');

% Define the initial mixture composition
mix1 = Mixture(system);
set(mix1, {'CH4'}, 'fuel', 1);
set(mix1, {'N2', 'O2'}, 'oxidizer', [79/21, 1]);

% Define the thermodynamic state
setTemperature(mix1, 300, 'K');
setPressure(mix1, 1, 'bar');
setEquivalenceRatio(mix1, 1); % set(mix1, 'equivalenceRatio', 0.5:0.01:2);

% Create the equilibrium solver
solver = EquilibriumSolver('flag_results', false);

% Solve chemical transformation at equilibrium
% tic
for i = N:-1:1
    mix2{i} = mix1.copy().setTemperature(temperature(i), 'K');
    solver.solve(mix2{i});
    % print_mixture(mix1, mix2{i});
end
% toc

% % feature('numcores');
% workers = 5;
% delete(gcp('nocreate'))
% parpool(workers);

% tic
% parfor i=1:1000
%     mix2{i} = mix1.copy().setTemperature(temperature(i), 'K');
%     solver.solve(mix2{i});
%     % print_mixture(mix1, mix2{i});
% end
% toc

% mix1 and mix2 are struct arrays of Mixture objects, each element of the array corresponds to a different equivalence ratio
% mix1 contains the initial mixture, mix2 contains the mixture at equilibrium
% Display the results
toc

plotFigure('T', cell2vector(mix2, 'T'), 'h', cell2vector(mix2, 'h'));
%% Same but using Array
% Example
% Import the necessary classes ( included in edit(fullfile(userpath,'startup.m')) )
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

clear

tic

% Definitions
temperature = linspace(1000, 6000, 1000);
N = length(temperature);

% Load NASA's database
DB = NasaDatabase();

% Define the chemical system, i.e., the set of species that may be present in the chemical transformation at equilibrium
system = ChemicalSystem(DB, 'Soot formation');

% Define the initial mixture composition
mix = Mixture(system);
set(mix, {'CH4'}, 'fuel', 1);
set(mix, {'N2', 'O2'}, 'oxidizer', [79/21, 1]);

% Define the thermodynamic state
mixArray = mix.setProperties('temperature', temperature, 'pressure', 1, 'equivalenceRatio', 1);

% Create the equilibrium solver
solver = EquilibriumSolver('flag_results', false);

% Solve chemical transformation at equilibrium
solver.solveArray(mixArray);

toc

plotFigure('T', [mixArray.T], 'h', [mixArray.T]);



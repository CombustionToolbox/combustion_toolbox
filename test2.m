% Example
% Import the necessary classes ( included in edit(fullfile(userpath,'startup.m')) )
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.EquilibriumSolver

tic
phi = linspace(0.5, 4, 1000);
nfrec = 10;

% Load NASA's database
DB = NasaDatabase();

% Define the chemical system, i.e., the set of species that may be present in the chemical transformation at equilibrium
system = ChemicalSystem(DB, 'Soot formation'); % HC/O2/N2

% Define the initial mixture composition
mix1 = Mixture(system);
set(mix1, {'CH4'}, 'fuel', 1);
set(mix1, {'N2', 'O2'}, 'oxidizer', [79/21, 1]);

% Define the thermodynamic state
setTemperature(mix1, 300, 'K');
setPressure(mix1, 1, 'bar');

% Create the equilibrium solver
solver = EquilibriumSolver('problemType', 'EV');

% Solve chemical transformation at equilibrium
% tic
numCases = length(phi);
mix1cell{numCases} = mix1.setEquivalenceRatio(phi(numCases));
mix2cell{numCases} = mix1cell{numCases}.copy();
solver.solve(mix2cell{numCases});
print(mix1cell{numCases}, mix2cell{numCases});

for i = numCases - 1:-1:1
    mix1cell{i} = mix1.copy().setEquivalenceRatio(phi(i));
    mix2cell{i} = mix1cell{i}.copy();
    solver.solve(mix2cell{i}, mix2cell{i + 1});
    if ~mod(i, nfrec) || i == 1
        print(mix1cell{i}, mix2cell{i});
    end

end
% toc

% % feature('numcores');
% workers = 5;
% delete(gcp('nocreate'))
% parpool(workers);

% tic
% parfor i=1:numCases
%     mix1cell{i} = mix1.setEquivalenceRatio(phi(i));
%     mix2cell{i} = mix1cell{i}.setTemperature(300, 'K');
%     mix2cell{i} = solver.solve(mix2cell{i});
%     % print(mix1cell{i}, mix2cell{i});
% end
% toc

% mix1 and mix2 are struct arrays of Mixture objects, each element of the array corresponds to a different equivalence ratio
% mix1 contains the initial mixture, mix2 contains the mixture at equilibrium
% Display the results

toc

plot_figure('phi', phi, 'T', cell2vector(mix2cell, 'T'));
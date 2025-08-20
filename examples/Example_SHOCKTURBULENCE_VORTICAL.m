% Example of a Vortical Shock-Turbulence Interaction

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.shockturbulence.*
import combustiontoolbox.utils.display.*

% Definitions
mach = combustiontoolbox.utils.clusteredMesh1D([1, 1.2], [1.2, 10], 32, 70); mach(1) = [];
FLAG_TCHEM_FROZEN = false;
FLAG_FROZEN = false;

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79/21, 1]);

% Define propertiescle
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'mach', mach);

% Get jump conditions
solver = JumpConditionsSolver('FLAG_TCHEM_FROZEN', FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', FLAG_FROZEN, 'FLAG_PAPER', false, 'FLAG_RESULTS', false);

% Solve problem
jumpData = solver.solve(mixArray);

% Invoke ShockTurbulenceSolver and select problem
shockTurbulence = ShockTurbulenceSolver('problemType', 'vortical', 'FLAG_TCHEM_FROZEN', FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', FLAG_FROZEN);

% Solve LIA
results = shockTurbulence.solve(mixArray, 'jumpconditions', jumpData);

% Report results
shockTurbulence.report(results);
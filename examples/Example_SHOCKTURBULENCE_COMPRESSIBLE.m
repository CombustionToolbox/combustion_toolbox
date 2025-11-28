% Example of an Acoustic Shock-Turbulence Interaction

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockturbulence.*
import combustiontoolbox.utils.display.*

% Definitions
mach = combustiontoolbox.utils.clusteredMesh1D([1, 1.2], [1.2, 10], 32, 70); mach(1) = [];
FLAG_TCHEM_FROZEN = true;
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

% Invoke ShockTurbulenceSolver and select problem
shockTurbulence = ShockTurbulenceSolver('problemType', 'compressible', 'FLAG_TCHEM_FROZEN', FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', FLAG_FROZEN);

% Solve LIA
results = shockTurbulence.solve(mixArray, 'eta', 0.1);

% Report results
shockTurbulence.report(results);
% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_LIMIT_RR_ALTITUDE
%
% Compute limit regular reflections at different altitudes considering a
% thermochemical frozen gas, a chemically frozen gas, and dissociation,
% ionization, vibrational excitation and electronic excitation. The
% calculations are carried out for a set of pre-shock Mach numbers 
% M1 = [1.75:0.1:20]
%    
% Air_ions == {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', ...
%              'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', ...
%              'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', ...
%              'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'C3O2', 'N', ...
%              'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', ...
%              'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
%              'N2minus', 'NCN', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
%              'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
%              'O2minus', 'O3'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% Note: This example uses the Standard Atmosphere Functions package [1]
%
% [1] Sky Sartorius (2024). Standard Atmosphere Functions. Available at
%     https://github.com/sky-s/standard-atmosphere, GitHub.
%     Retrieved June 12, 2024.
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Jun 11 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*
import combustiontoolbox.common.Units

% Definitions
S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
N_Oxidizer = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
machNumber = 1.75:0.1:20.5;
z = [0, 15000, 30000];
numCases = length(machNumber);
numAltitude = length(z);

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air_ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2', 'Ar', 'CO2'}, [78.084, 20.9476, 0.9365, 0.0319] / 20.9476);

% Initialization of Mixture objects
mixArray1 = repmat(mix, [numCases, 3, numAltitude]);
mixArray2 = repmat(mix, [numCases, 3, numAltitude]);
mixArray2_1 = repmat(mix, [numCases, 3, numAltitude]);
mixArray3 = repmat(mix, [numCases, 3, numAltitude]);

% Initialization of EquilibriumSolver
equilibriumSolver = EquilibriumSolver();

for i = 3:3

    switch i
        case 1
            equilibriumSolver.FLAG_FROZEN = false;
            equilibriumSolver.FLAG_TCHEM_FROZEN = false;
        case 2
            equilibriumSolver.FLAG_FROZEN = true;
            equilibriumSolver.FLAG_TCHEM_FROZEN = false;
        case 3
            equilibriumSolver.FLAG_FROZEN = false;
            equilibriumSolver.FLAG_TCHEM_FROZEN = true;
    end

    for j = numAltitude
        % Get temperature and pressure
        [~, ~, temperature, pressure] = atmos(z(j));
        pressure = Units.convert(pressure, 'Pa', 'bar');
        
        % Define properties
        mixArray1(:, i, j) = setProperties(mix, 'temperature', temperature, 'pressure', pressure, 'M1', machNumber);

        % Initialize solver
        solver = ShockSolver('problemType', 'SHOCK_POLAR_LIMITRR', 'equilibriumSolver', equilibriumSolver, 'tol0', 1e-7);
        
        % Solve problem
        [mixArray1(:, i, j),...
         mixArray2(:, i, j),...
         mixArray2_1(:, i, j),...
         mixArray3(:, i, j)] = solver.solveArray(mixArray1(:, i, j));
    end

end

% Definitions
ax = setFigure();
colors = [175 221 233; 95 188 211; 0 102 128] / 255;

for j = 1:numAltitude
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Wave angle limit RR [deg]', [mixArray2_1(:, i, j).beta], 'linestyle', 'k-', 'ax', ax, 'color', colors(j, :));
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Wave angle limit RR [deg]', [mixArray2_1(:, i, j).beta], 'linestyle', 'k--', 'ax', ax, 'color', colors(j, :));
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Wave angle limit RR [deg]', [mixArray2_1(:, i, j).beta], 'linestyle', 'k:', 'ax', ax, 'color', colors(j, :));
end

for j = 1:numAltitude
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Deflection angle limit RR [deg]', [mixArray2_1(:, i, j).theta], 'linestyle', 'k-', 'ax', ax, 'color', colors(j, :));
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Deflection angle limit RR [deg]', [mixArray2_1(:, i, j).theta], 'linestyle', 'k--', 'ax', ax, 'color', colors(j, :));
    ax = plotFigure('Pre-shock Mach number', [mixArray1(:, i, j).mach], 'Deflection angle limit RR [deg]', [mixArray2_1(:, i, j).theta], 'linestyle', 'k:', 'ax', ax, 'color', colors(j, :));
end
% -------------------------------------------------------------------------
% EXAMPLE: TP_scoggins2015
%
% Compute equilibrium composition for a defined set of temperatures
% (200, 5000) [K] at atmospheric pressure (1.01325 bar) for a Si-C6H5OH_phenol
% mixture at standard conditions, and considering a set of 177 species 
% (40 in condensed phase)
% 
% This example is obtained from [1]
% 
% [1] Scoggins, J. B., & Magin, T. E. (2015). Gibbs function continuation
%     for linearly constrained multiphase equilibria. Combustion and Flame,
%     162(12), 4514-4522.
%
% @author: Alberto Cuadra Lara
%                 
% Last update Oct 06 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*

% Definitions
plotConfig = PlotConfig('mintolDisplay', 1e-2);

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'Si', 'C6H5OH_phenol'}, [1, 9]);

% Define properties
mixArray = setProperties(mix, 'temperature', 200:10:5000, 'pressure', 1 * 1.01325);

% Initialize solver
solver = EquilibriumSolver('problemType', 'TP', 'plotConfig', plotConfig);

% Solve problem
solver.solveArray(mixArray);

% Generate report
report(solver, mixArray);
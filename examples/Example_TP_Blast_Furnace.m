% -------------------------------------------------------------------------
% EXAMPLE: TP_furnance
%
% Compute equilibrium composition of a blast furnace at defined temperature
% 1050 K and pressure 1.01325 bar. A set of 10 gaseous species and 2 
% condensed species are considered
%
% This example is obtained from [1]
% 
% [1] W.D. Madeley & J.M. Toguri (1973) The application of free energy
%     minimization techniques to determine equilibrium compositions in
%     systems of metallurgical interest, Canadian Metallurgical Quarterly,
%     12:1, 71-78. DOI: 10.1179/cmq.1973.12.1.71
%   
% Species == {'O2', 'N2', 'H2O', 'CH4', 'CO', 'CO2', 'H2', ...
%             'CHO_M', 'CH2O_M', 'OH', 'Febab', 'CaObcrb'}
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Jun 04 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.findIndex

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'O2', 'N2', 'H2O', 'CH4', 'Fe3O4bcrb', 'Febab', 'Cbgrb', 'CaCO3bcrb', 'CaObcrb'},...
    [20.46, 187.1, 1.775, 2.2554, 13.1, 3.527, 85.59, 0.1499, 0.6063]);

% Define properties
mixArray = setProperties(mix, 'temperature', 1050, 'pressure', 1 * 1.01325);

% Initialize solver
solver = EquilibriumSolver('problemType', 'TP', 'tolMoles', 1e-25, 'FLAG_RESULTS', true);

% Solve problem
solver.solveArray(mixArray);

% Print results
displaySpecies = {'O2', 'N2', 'H2O', 'CH4', 'CO', 'CO2', 'H2', 'OH', 'Febab', 'CaObcrb'};
indexSpecies = findIndex(system.listSpecies, displaySpecies);
moles_RAND = [1.793e-19; 1.871e2; 2.377e-1; 1.97e-3; 8.143e1; 6.865; 6.641; 8.414e-11; 4.283e1; 7.562e-1];
moles_CT = mixArray.Xi(indexSpecies) * mixArray.N;

% Display table with species and molar coposition
T = table(system.listSpecies(indexSpecies)', moles_CT, moles_RAND, moles_RAND - moles_CT);
T.Properties.VariableNames = {'Species', 'Moles', 'Moles RAND', 'Error [moles]'}


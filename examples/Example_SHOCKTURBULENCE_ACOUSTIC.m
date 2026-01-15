% -------------------------------------------------------------------------
% EXAMPLE: SHOCKTURBULENCE_ACOUSTIC
%
% This example computes the Linear Interaction Analysis (LIA) solution for
% a shock–turbulence interaction involving upstream acoustic fluctuations
% only, with no vortical or entropic disturbances present in the incoming
% flow.
%
% The upstream turbulent field is modeled as a purely dilatational acoustic
% mode, consisting of propagating pressure and velocity fluctuations that
% satisfy the linearized isentropic relations. No solenoidal (vortical)
% component or correlated entropic density fluctuations are prescribed.
%
% The interaction with a normal shock is analyzed under the assumptions of
% linear perturbations, inviscid flow, and thermochemical equilibrium across
% a thin relaxation layer.
%
% The formulation follows the LIA framework described in:
%
%   Cuadra, A., Williams, C. T., Di Renzo, M., & Huete, C., The role of
%   compressibility and vibrational-excitation in hypersonic shock–turbulence
%   interactions, Journal of Fluid Mechanics (under review).
%
% In this acoustic configuration, all upstream dilatational content is
% associated with propagating acoustic waves. Consequently, the entropic–
% vortical correlation parameter :math:`\chi` is not defined, and the
% vortical turbulent kinetic energy is identically zero.
%
% @author: Alberto Cuadra Lara
%
% Last update: December 16, 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockturbulence.*
import combustiontoolbox.utils.display.*

% Definitions
mach = combustiontoolbox.utils.clusteredMesh1D([1, 5], [5, 10], 32, 70); mach(1) = [];

% Define caloric gas model
caloricGasModel = CaloricGasModel.imperfect;

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
shockTurbulence = ShockTurbulenceSolver('problemType', 'acoustic', 'caloricGasModel', caloricGasModel);

% Solve LIA
[averages, mixArray1, mixArray2] = shockTurbulence.solve(mixArray);

% Report results
shockTurbulence.report(averages, mixArray1, mixArray2);
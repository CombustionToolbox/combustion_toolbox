% -------------------------------------------------------------------------
% Example of a vortical Shock–Turbulence Interaction using LIA
%
% This example computes the Linear Interaction Analysis (LIA) solution for
% a shock–turbulence interaction involving upstream vortical fluctuations
% only, with no entropic or acoustic disturbances present in the incoming
% flow.
%
% The upstream turbulent field is modeled as a purely solenoidal vortical
% mode, characterized by velocity fluctuations with zero dilatation and
% no associated density or pressure perturbations. As a result, the
% incoming flow is incompressible in the LIA sense.
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
% In this vortical configuration, the upstream turbulent kinetic energy is
% entirely solenoidal, such that
%
% .. math::
%
%    \mathrm{TKE}_1 = \mathrm{TKE}_{1,r}
%
% and no dilatational contribution is present. Consequently, the
% dilatational-to-solenoidal TKE ratio satisfies
%
% .. math::
%
%    \eta = 0
%
% and the entropic–vortical correlation parameter :math:`\chi` is not
% defined.
%
% @author: Alberto Cuadra Lara
%
% Last update: December 16, 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.shockturbulence.*
import combustiontoolbox.utils.display.*

% Definitions
mach = combustiontoolbox.utils.clusteredMesh1D([1, 1.2], [1.2, 10], 32, 70); mach(1) = [];

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
mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'mach', mach);

% Invoke ShockTurbulenceSolver and select problem
shockTurbulence = ShockTurbulenceSolver('problemType', 'vortical', 'caloricGasModel', caloricGasModel);

% Solve LIA
[averages, mixArray1, mixArray2] = shockTurbulence.solve(mixArray1);

% Report results
shockTurbulence.report(averages, mixArray1, mixArray2);
% -------------------------------------------------------------------------
% Example of a vortical-entropic Shock-Turbulence Interaction using LIA
%
% This example computes the Linear Interaction Analysis (LIA) solution for
% a shock–turbulence interaction involving upstream vortical–entropic
% fluctuations only, with no propagating acoustic content.
%
% The upstream turbulent field is modeled as a solenoidal vortical mode with
% correlated entropic (dilatational) density fluctuations. In this formulation,
% all dilatational content is implicit in the vortical–entropic disturbance and
% no independent acoustic (traveling-wave) fluctuations are prescribed.
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
% The coupling between entropic (dilatational) density fluctuations and
% vortical velocity perturbations is prescribed through the correlation
% parameter :math:`\chi`, defined according to Eq. (3.9) of the reference 
% paper as
%
% .. math::
%
%    \langle \chi \rangle =
%    \frac{\langle \delta \rho_1^{e} \, \delta u_1^{r} \rangle}
%         {\langle (\delta u_1^{r})^2 \rangle}
%    \frac{\langle c_1 \rangle}{\langle \rho_1 \rangle}
%
% where :math:`\delta \rho_1^{e}` denotes entropic density fluctuations and
% :math:`\delta u_1^{r}` denotes rotational velocity fluctuations ahead of
% the shock.
%
% In this example, :math:`\chi = -0.1`, corresponding to negatively correlated
% entropic and vortical disturbances. This choice represents a compressible
% vortical field in which density and velocity fluctuations are anti-correlated,
% enhancing baroclinic vorticity generation at the shock.
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
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'mach', mach);

% Invoke ShockTurbulenceSolver and select problem
shockTurbulence = ShockTurbulenceSolver('problemType', 'vortical_entropic', 'caloricGasModel', caloricGasModel);

% Solve LIA
results = shockTurbulence.solve(mixArray, 'chi', -0.1);

% Report results
shockTurbulence.report(results);
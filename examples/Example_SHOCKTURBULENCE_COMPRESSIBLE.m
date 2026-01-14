% -------------------------------------------------------------------------
% Example of compressible Shock-Turbulence Interaction using LIA
%
% This example computes the Linear Interaction Analysis (LIA) solution for 
% a compressible shock–turbulence interaction involving upstream vortical
% and acoustic fluctuations.
% 
% The upstream turbulent field is modeled as a superposition of solenoidal
% (vortical) and dilatational modes. The dilatational content is partitioned
% into:
%
%   (i) entropic fluctuations correlated with vortical disturbances, and
%  (ii) acoustic (traveling-wave) fluctuations.
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
% The relative contribution of acoustic (dilatational) fluctuations to the
% upstream turbulent kinetic energy is prescribed through the parameter
% :math:`\eta`, defined as the ratio of dilatational to solenoidal TKE,
% 
% .. math::
% 
%    \eta = \frac{\mathrm{TKE}_{1,a}}{\mathrm{TKE}_{1,r}}
% 
% where subscripts :math:`a` and :math:`r` denote acoustic and rotational
% (vortical–entropic) components, respectively.
% 
% In this example, the compressible case corresponds to
% 
% .. math::
% 
%    \eta = 0.1
% 
% indicating that 10% of the upstream turbulent kinetic energy is associated
% with acoustic (dilatational) fluctuations.
%
% The correlation parameter :math:`\chi` characterizes the entropic (dilatational)
% density fluctuations that are correlated with vortical disturbances.
% Although chi also represents dilatational content, it is implicit in the
% vortical–entropic mode and does not correspond to propagating acoustic
% energy.
%
% In this example, :math:`\chi` is set to zero, corresponding to vortical
% fluctuations without correlated entropic disturbances.
%
% @author: Alberto Cuadra Lara
%                 
% Last update January 12, 2026
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockturbulence.*
import combustiontoolbox.utils.display.*

% Definitions
mach = combustiontoolbox.utils.clusteredMesh1D([1, 1.2], [1.2, 10], 32, 70); mach(1) = [];
caloricGasModel = CaloricGasModel.imperfect;

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79/21, 1]);

% Define properties
mixArray = setProperties(mix, 'temperature', 300, 'pressure', 1 * 1.01325, 'mach', mach, 'eta', 0.1);

% Invoke ShockTurbulenceSolver and select problem
shockTurbulence = ShockTurbulenceSolver('problemType', 'compressible', 'caloricGasModel', caloricGasModel);

% Update viscosity model
shockTurbulence.shockTurbulenceModel.viscosityModel = 'sutherland';

% Solve LIA
[averages, mixArray1, mixArray2] = shockTurbulence.solve(mixArray);

% Report results
shockTurbulence.report(averages, mixArray1, mixArray2);
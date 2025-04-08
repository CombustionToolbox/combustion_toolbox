% -------------------------------------------------------------------------
% EXAMPLE: Helmholtz decomposition
%
% Perform Helmholtz decomposition of a velocity field and analyze the
% turbulence spectra of the original, solenoidal, and dilatational fields.
%
% Note: The data used is a downsampled velocity field from a Direct Numerical
% Simulation (DNS) of a Homogeneous Isotropic Turbulence (HIT) case carried
% out with the HTR solver [1]. Take the data as an example to illustrate
% the use of the HelmholtzSolver and TurbulenceSpectra classes.
%
% References:
% [1] Di Renzo, M., Fu, L., & Urzay, J. (2020). HTR solver: An
%     open-source exascale-oriented task-based multi-GPU high-order
%     code for hypersonic aerothermodynamics. Computer Physics
%     Communications, 255, 107262.
%
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Nov 23 2024
% -------------------------------------------------------------------------

% Import required packages
import combustiontoolbox.turbulence.*

% Load field variables
filename = 'velocityField.h5';
filePath = which(filename);
rho = double(h5read(filePath, ['/', 'density']));
u = double(h5read(filePath, ['/', 'velocity/u']));
v = double(h5read(filePath, ['/', 'velocity/v']));
w = double(h5read(filePath, ['/', 'velocity/w']));

% Convert velocity to VelocityField object
velocity = VelocityField(u, v, w);

% Initialize HelmholtzSolver
solver = HelmholtzSolver();

% Perform Helmholtz decomposition
[solenoidal, dilatational, STOP] = solver.solve(velocity, 'density', rho);

% Compute turbulent kinetic energy (TKE)
K_solenoidal = solenoidal.getTurbulentKineticEnergy(rho);
K_dilatational = dilatational.getTurbulentKineticEnergy(rho);

% Get dilatational over solenoidal TKE ratio
eta = K_dilatational / K_solenoidal;

% Analyze turbulence spectra using TurbulenceSpectra
analyzer = TurbulenceSpectra();

% Compute energy spectra for the original, solenoidal, and dilatational fields
[EK1, k1] = analyzer.getEnergySpectra(velocity);     % Original field
[EK2, k2] = analyzer.getEnergySpectra(solenoidal);   % Solenoidal component
[EK3, k3] = analyzer.getEnergySpectra(dilatational); % Dilatational component

% Plot results
analyzer.plot(k1, EK1, EK2, EK3);
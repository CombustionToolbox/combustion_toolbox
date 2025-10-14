% -------------------------------------------------------------------------
% COMBUSTION TOOLBOX @v1.2.7
% A MATLAB-GUI based open-source tool for solving gaseous combustion problems
%
% Type of problems:
%   * TP ------------------> Equilibrium composition at defined T and p
%   * HP ------------------> Adiabatic T and composition at constant p
%   * SP ------------------> Isentropic compression/expansion to a specified p
%   * TV ------------------> Equilibrium composition at defined T and constant v
%   * EV ------------------> Adiabatic T and composition at constant v
%   * SV ------------------> Isentropic compression/expansion to a specified v
%   * SHOCK_I -------------> Planar incident shock wave
%   * SHOCK_R -------------> Planar reflected shock wave
%   * SHOCK_OBLIQUE -------> Oblique incident shock wave
%   * SHOCK_OBLIQUE_R -----> Oblique incident and reflected states
%   * SHOCK_POLAR ---------> Shock polar diagrams
%   * SHOCK_POLAR_R -------> Shock polar diagrams for incident and reflected states
%   * SHOCK_POLAR_LIMITRR -> Shock polar diagrams in the limit of regular reflection
%   * SHOCK_IDEAL_GAS -----> Planar incident shock wave for a fixed adiabatic index
%   * DET -----------------> Chapman-Jouguet Detonation
%   * DET_R ---------------> Reflected Chapman-Jouguet Detonation
%   * DET_OBLIQUE ---------> Oblique Detonation
%   * DET_POLAR -----------> Detonation polar diagrams
%   * DET_OVERDRIVEN ------> Over-driven Detonation    
%   * DET_OVERDRIVEN_R ----> Over-driven reflected Detonation
%   * DET_UNDERDRIVEN -----> Under-driven Detonation
%   * DET_UNDERDRIVEN_R ---> Under-driven reflected Detonation
%   * ROCKET --------------> Propellant rocket performance
%   * HELMHOLTZ -----------> Helmholtz-Hodge decomposition of a velocity field
%   * SPECTRA -------------> Turbulence spectra analysis
%
% SEE THE EXAMPLES OR WEBSITE TO KNOW HOW TO START USING THE COMBUSTION TOOLBOX
%
% WEBSITE: https://combustion-toolbox-website.readthedocs.io/ 
%          or type in the promt "websiteCT".
%
% Please to send feedback or inquiries type in the promt "uifeedback".
%
% Thank you for using the Combustion Toolbox!
%
% Citing:
%     Cuadra, A., Huete, C., & Vera, M. (2024). Combustion Toolbox: An
%     open-source thermochemical code for gas- and condensed-phase
%     problems involving chemical equilibrium. arXiv:2409.15086.
%
%     Cuadra, A., Huete, C., & Vera, M. (2025). Combustion Toolbox: A
%     MATLAB-GUI based open-source tool for solving gaseous combustion
%     problems. (v1.2.7). Zenodo. https://doi.org/10.5281/zenodo.5554911.
%
% @author: Alberto Cuadra Lara
%                  
% Last update October 14 2025
% -------------------------------------------------------------------------
help CONTENTS.m

% Set path
INSTALL();
% Display splash
gui_display_splash('pause', 2);
% Check for updates
combustiontoolbox.utils.checkUpdate();
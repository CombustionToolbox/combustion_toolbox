% -------------------------------------------------------------------------
% COMBUSTION TOOLBOX @v0.9.0
% A MATLAB-GUI based open-source tool for solving gaseous combustion problems
%
% Type of problems:
%   * TP -----------------> Equilibrium composition at defined T and p
%   * HP -----------------> Adiabatic T and composition at constant p
%   * SP -----------------> Isentropic compression/expansion to a specified p
%   * TV -----------------> Equilibrium composition at defined T and constant v
%   * EV -----------------> Adiabatic T and composition at constant v
%   * SV -----------------> Isentropic compression/expansion to a specified v
%   * SHOCK_I ------------> Planar incident shock wave
%   * SHOCK_R ------------> Planar reflected shock wave
%   * SHOCK_OBLIQUE ------> Oblique incident shock wave
%   * SHOCK_OBLIQUE_R ----> Oblique incident and reflected states
%   * SHOCK_POLAR --------> Shock polar plots
%   * SHOCK_POLAR_R ------> Shock polar plots for incident and reflected states
%   * SHOCK_IDEAL_GAS ----> Planar incident shock wave for a fixed adibatic index
%   * DET ----------------> Chapman-Jouguet Detonation
%   * DET_R --------------> Reflected Chapman-Jouguet Detonation
%   * DET_OBLIQUE --------> Oblique Detonation
%   * DET_OVERDRIVEN -----> Overdriven Detonation    
%   * DET_OVERDRIVEN_R ---> Overdriven reflected Detonation    
%   * ROCKET -------------> Propellant rocket performance   
%
% SEE THE EXAMPLES OR WIKI TO KNOW HOW TO START USING COMBUSTION TOOLBOX 
% LIST OF TUTORIAL SCRIPTS:
%   * Example_TP
%   * Example_HP
%   * Example_HP_PROPELLANTS
%   * Example_HP_MIXTEMP
%   * Example_SP
%   * Example_TV
%   * Example_EV
%   * Example_SV
%   * Example_SV_FROZEN
%   * Example_SHOCK_I
%   * Example_SHOCK_I_IONIZATION
%   * Example_SHOCK_R
%   * Example_SHOCK_OBLIQUE_BETA
%   * Example_SHOCK_OBLIQUE_THETA
%   * Example_SHOCK_OBLIQUE_R
%   * Example_SHOCK_POLAR
%   * Example_SHOCK_POLAR_R
%   * Example_DET
%   * Example_DET_R
%   * Example_DET_OVERDRIVEN
%   * Example_DET_OVERDRIVEN_R
%   * Example_ROCKET
%
% Please to send feedback or inquiries run uifeedback
% Thank you for testing Combustion Toolbox!
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update March 24 2022
% -------------------------------------------------------------------------
help Combustion_Toolbox.m

% INDICATE FILES ON PATH
Combustion_Toolbox_setPath;
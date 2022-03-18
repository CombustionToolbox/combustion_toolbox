% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_OBLIQUE
%
% Compute pre-shock and post-shock state for a oblique incident shock wave
% at standard conditions, a set of 20 species considered and a initial 
% shock front velocities u1 = 1500 [m/s]
%    
% Air == {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4',...
%         'N3','C','CO','CO2','Ar','H2O','H2','H','He'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update March 18 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air_ions');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2'};
self.PD.proportion_inerts_O2 = 79/21;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
overdriven = 2:2:14;
self = set_prop(self, 'u1', 3.477205457292110e+02 * overdriven);
%% SOLVE PROBLEM
self = SolveProblem(self, 'SHOCK_OBLIQUE');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
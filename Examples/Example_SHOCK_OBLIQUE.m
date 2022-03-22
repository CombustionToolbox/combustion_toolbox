% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_OBLIQUE
%
% Compute pre-shock and post-shock state for a oblique incident shock wave
% at standard conditions, a set of 20 species considered, a initial 
% shock front velocities u1 = a1 * 10 [m/s], and a deflection angle 
% theta = 20 [deg]
%    
% Air_ions == {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3',...
%              'N2O4','N3','eminus','Nminus','Nplus','NOplus','NO2minus',...
%              'NO3minus','N2plus','N2minus','N2Oplus','Oplus','Ominus',...
%              'O2plus', 'O2minus,'CO2','CO','COplus','C','Cplus',...
%              'Cminus','CN','CNplus','CNminus','CNN','NCO','NCN','Ar',...
%              'Arplus'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update March 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air_ions');
% self = App({'O2', 'N2', 'Ar', 'CO2'}); % Frozen
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
overdriven = 10;
self = set_prop(self, 'u1', 3.472107491008314e+02 * overdriven, 'theta', 20);
%% SOLVE PROBLEM
self = SolveProblem(self, 'SHOCK_OBLIQUE');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
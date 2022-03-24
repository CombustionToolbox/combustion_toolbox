% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED
%
% Compute shock polar plots at standard conditions, a set of 39 species
% considered, and a set of initial shock front velocities u1/a1 = [2:2:14]
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
% Last update March 21 2021
% -------------------------------------------------------------------------

%% INITIALIZE
% self = App('Air_ions');   
self = App({'O2', 'N2', 'Ar', 'CO2'}); % Frozen
% self = App({'O2'}); % Frozen
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
% range1 = logspace(0, 1, 300); range1 = range1(range1 < 5);
% overdriven = [range1, linspace(5, 14, 30)]; overdriven = overdriven(overdriven > 1);
% overdriven = 2:2:14;
% overdriven = [2, 3, 5, 14];
overdriven = 5;
self = set_prop(self, 'u1', 3.472107491008314e+02 * overdriven, 'theta', 4);
%% SOLVE PROBLEM
self = SolveProblem(self, 'SHOCK_POLAR_R');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
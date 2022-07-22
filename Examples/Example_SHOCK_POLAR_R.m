% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_REFLECTED
%
% Compute shock polar plots at standard conditions, a set of 39 species
% considered, an initial shock front Mach numbers = 6, and a deflection 
% angle theta = 25;
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
% Last update July 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air_ions');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.N_Oxidizer = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
Mach_number = 6;
self = set_prop(self, 'u1', 301.8203 * Mach_number, 'theta', 25);
%% SOLVE PROBLEM
self = SolveProblem(self, 'SHOCK_POLAR_R');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
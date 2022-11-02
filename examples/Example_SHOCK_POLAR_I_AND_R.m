% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_INCIDENT_AND_REFLECTED
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
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update July 22 2022
% -------------------------------------------------------------------------

%% SOLVE REFLECTED SHOCK
% INITIALIZE
self = App({'N2', 'O2', 'Ar'});
% INITIAL CONDITIONS
self = set_prop(self, 'TR', 226.51, 'pR', 0.01197);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar'};
self.PD.N_Oxidizer = [78, 21, 1] ./ 21;
% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
Mach_number = 20;
self = set_prop(self, 'u1', 301.8203 * Mach_number, 'theta', 35);
% SOLVE PROBLEM
self = solve_problem(self, 'SHOCK_POLAR_R');
% DISPLAY RESULTS (PLOTS)
post_results(self);
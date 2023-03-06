% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR
%
% Compute shock polar plots at standard conditions, a set of 39 species
% considered, and a set of initial shock front velocities u1/a1 = [2, 3, 5, 14]
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
% Last update Jan 10 2023
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air_ions');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.N_Oxidizer = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
sound_velocity = compute_sound(self.PD.TR.value, self.PD.pR.value, self.PD.S_Oxidizer, self.PD.N_Oxidizer, 'self', self);
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
Mach_number = [2, 3, 5, 14];
self = set_prop(self, 'u1', sound_velocity * Mach_number);
%% SOLVE PROBLEM
self = solve_problem(self, 'SHOCK_POLAR');
%% DISPLAY RESULTS (PLOTS)
post_results(self);
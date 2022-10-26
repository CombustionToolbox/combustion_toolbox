% -------------------------------------------------------------------------
% EXAMPLE: TP_class_ex9b
%
% Compute equilibrium composition for a defined temperature (3000 K) and a
% set of pressures (1.01325 to 101.325 bar) for an air mixture at standard 
% conditions, and a set of 4 species is considered.
%   
% LS == {'N2', 'O2', 'O', 'NO'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update July 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'N2', 'O2', 'O', 'NO'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'N2', 'O2'};
self.PD.N_Oxidizer = [0.79, 0.21];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', [1:1:100] * 1.01325, 'TP', 3000); 
%% SOLVE PROBLEM
self = solve_problem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
post_results(self);
%% PRINT MOLES
species = self.S.LS';
i = find(self.PD.pP.value == 10.1325);
molar_fraction = self.PS.strP{i}.Xi;
moles = molar_fraction * self.PS.strP{i}.N;
table(species, moles, molar_fraction)
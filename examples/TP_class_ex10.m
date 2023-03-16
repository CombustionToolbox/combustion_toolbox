% -------------------------------------------------------------------------
% EXAMPLE: TP_class_ex10
%
% Compute equilibrium composition for a defined temperature (2000 K) and 
% pressures (1.01325 bar) for a rich CH4-air mixture at
% standard conditions, a set of 6 species considered, and a equivalence
% ratio phi (1.25) [-]
%   
% LS_rich_complete == {'CO2', 'H2O', 'N2', 'CO', 'H2'}
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
self = App({'CO2', 'H2O', 'N2', 'CO', 'H2'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298, 'pR', 1 * 1.01325, 'phi', 1.25);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'N2', 'O2'};
self.PD.ratio_oxidizers_O2 = [79, 21] ./ 21;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', 1 * 1.01325, 'TP', 2000); 
%% SOLVE PROBLEM
self = solve_problem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
post_results(self);
%% PRINT MOLES
species = self.S.LS';
i = 1;
molar_fraction = self.PS.strP{i}.Xi;
moles = molar_fraction * self.PS.strP{i}.N;
table(species, moles, molar_fraction)
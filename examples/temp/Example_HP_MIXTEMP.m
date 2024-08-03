% -------------------------------------------------------------------------
% EXAMPLE: HP MIXTEMP
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich CH4-air mixtures at
% standard conditions except for the air which is at 380 K. Also, a set
% of 24 species considered and a set of equivalence ratios phi contained
% in (0.5, 5) [-]
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;

self.PD.T_Fuel     = 300;
self.PD.T_Oxidizer = [380, 380, 380, 380];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value); 
%% SOLVE PROBLEM
self = solve_problem(self, 'HP');
%% DISPLAY RESULTS (PLOTS)
post_results(self);
% -------------------------------------------------------------------------
% EXAMPLE: HP MIXTEMP
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich CH4-air mixtures at
% standard conditions except for the air which is at 380 K. Also, a set
% of 24 species considered and a set of equivalence ratios phi contained
% in (0.5, 5) [-]
%   
% Soot formation == {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
%                    'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
%                    'NO','HCO','NH2','NH','N','CH','Cbgrb'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Feb 19 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;

self.PD.T_Fuel     = 300;
self.PD.T_Oxidizer = 380;
self.PD.T_Inert    = [380, 380, 380];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'HP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
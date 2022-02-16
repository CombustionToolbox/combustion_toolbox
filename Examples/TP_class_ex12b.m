% -------------------------------------------------------------------------
% EXAMPLE: TP_class_ex12b
%
% Compute adiabatic temperature and equilibrium composition at constant
% volume for lean to rich CH4-air mixtures at standard conditions, a set
% of 24 species considered and a set of equivalence ratios phi contained
% in (0.7, 1.5) [-]
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
% Last update Feb 16 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298, 'pR', 1 * 1.01325, 'phi', 0.7:0.01:1.5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2'};
self.PD.proportion_inerts_O2 = 79/21;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'EV');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
%% PRINT MOLES
species = self.S.LS';
i = find(self.PD.phi.value == 1.1);
molar_fraction = self.PS.strP{i}.Xi;
moles = molar_fraction * self.PS.strP{i}.N;
table(species, moles, molar_fraction)
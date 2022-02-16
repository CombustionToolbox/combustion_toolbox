% -------------------------------------------------------------------------
% EXAMPLE: TP_class_ex12a
%
% Compute adiabatic temperature and equilibrium composition at constant
% volume for lean to rich CH4-air mixtures at standard conditions, assuming
% a complete combustion, and a set of equivalence ratios phi contained
% in (0.7, 1.5) [-]
%   
% LS == 'Complete'
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
self = App('Complete');
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
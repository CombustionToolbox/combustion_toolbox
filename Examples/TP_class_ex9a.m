% -------------------------------------------------------------------------
% EXAMPLE: TP_class_ex9a
%
% Compute equilibrium composition for a set of defined temperatures 
% (400, 5000 K) and pressure (1.01325 bar) for an air mixture at standard 
% conditions, and a set of 4 species is considered.
%   
% LS == {'N2', 'O2', 'O', 'NO'}
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
self = App({'N2', 'O2', 'O', 'NO'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 298, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2'};
self.PD.proportion_inerts_O2 = 79/21;
self.PD.N_Oxidizer = 0.21;
self.PD.N_Inert    = 0.79;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 400:100:5000); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
%% PRINT MOLES
species = self.S.LS';
i = find(self.PD.TP.value == 3000);
molar_fraction = self.PS.strP{i}.Xi;
moles = molar_fraction * self.PS.strP{i}.N;
table(species, moles, molar_fraction)
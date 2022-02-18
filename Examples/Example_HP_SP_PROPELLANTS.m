% -------------------------------------------------------------------------
% EXAMPLE: HP and SP PROPELLANTS
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (0.5, 5) [-]
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Oct 22 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('HYDROGEN_L');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 2);
self.PD.S_Fuel     = {'H2bLb'};
self.PD.S_Oxidizer = {'O2bLb'};
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'HP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
% -------------------------------------------------------------------------
%% INITIALIZE
self2 = App('HYDROGEN_L');
%% INITIAL CONDITIONS
self2 = set_prop(self2, 'TR', self.PS.strP{1}.T, 'pR', self.PS.strP{1}.p);
self2.PD.S_Oxidizer = self.S.LS;
self2.PD.N_Oxidizer = moles(self.PS.strP{1})';
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self2 = set_prop(self2, 'pP', 0.58059); 
%% SOLVE PROBLEM
self2 = SolveProblem(self2, 'SP');
%% DISPLAY RESULTS (PLOTS)
postResults(self2);
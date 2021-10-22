% -------------------------------------------------------------------------
% EXAMPLE: DET_OVERDRIVEN
%
% Compute pre-shock and post-shock state for a planar overdriven detonation
% considering Chapman-Jouguet (CJ) theory for a stoichiometric CH4-air
% mixture at standard conditions, a set of 24 species considered and a set
% of overdrives contained in (1,10) [-].
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
% Last update Oct 22 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot Formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 0.5:0.01:5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
overdriven = 1:0.1:10;
self = set_prop(self, 'overdriven', overdriven, 'phi', self.PD.phi.value(1) * ones(1, length(overdriven)));
% condition
%% SOLVE PROBLEM
self = SolveProblem(self, 'DET_OVERDRIVEN');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
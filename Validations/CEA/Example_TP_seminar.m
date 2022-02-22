% -------------------------------------------------------------------------
% EXAMPLE: TP_seminar
%
% Compute equilibrium composition at defined temperature (e.g., 1500 K) and
% pressure (e.g., 1.01325 bar) for a lean H2-Br2 mixture at
% standard conditions, and a set of 3 species considered.
%   
% species == {'H2', 'Br2', 'HBr'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Dec 10 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'H2', 'Br2', 'HBr'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Fuel     = {'H2', 'Br2'};
self.PD.N_Fuel     = [0.75, 0.25];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 1500); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
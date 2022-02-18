% -------------------------------------------------------------------------
% EXAMPLE: SP
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of pressure (1.01325, 1013.25 bar) for a rich CH4-air mixture
% at standard conditions, a set of 24 species considered, and a equivalence
% ratio phi 1.5 [-]
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
% Last update Feb 18 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1.5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', 1.01325 * logspace(0, 3, 200)); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'SP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
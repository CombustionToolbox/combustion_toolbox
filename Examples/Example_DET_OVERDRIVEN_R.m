% -------------------------------------------------------------------------
% EXAMPLE: DET OVERDRIVEN REFLECTED
%
% Compute pre-shock and post-shock state for a reflected planar overdriven
% detonation considering Chapman-Jouguet (CJ) theory for a stoichiometric 
% CH4-air mixture at standard conditions, a set of 26 species considered
% and a set of overdrives contained in (1,10) [-].
%   
% Soot formation == {'CO2','CO','H2O','H2','O2','N2','He','Ar','Cbgrb',...
%                    'C2','C2H4','CH','CH','CH3','CH4','CN','H',...
%                    'HCN','HCO','N','NH','NH2','NH3','NO','O','OH'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 22 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot Formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.ratio_oxidizers_O2 = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
overdriven = 1:0.1:10;
self = set_prop(self, 'overdriven', overdriven);
%% SOLVE PROBLEM
self = SolveProblem(self, 'DET_OVERDRIVEN_R');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
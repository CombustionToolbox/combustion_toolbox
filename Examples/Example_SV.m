% -------------------------------------------------------------------------
% EXAMPLE: SV
% Compute Isentropic compression/expansion and equilibrium composition at 
% a defined set of volume ratios (0.5, 2) for a lean CH4-air mixture at
% 700 K and 10 bar, a set of 26 species considered, and a equivalence
% ratio phi 0.5 [-]
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
% Last update Feb 19 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Soot formation');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 700, 'pR', 10, 'phi', 0.5);
self.PD.S_Fuel     = {'CH4'};
self.PD.S_Oxidizer = {'O2'};
self.PD.S_Inert    = {'N2', 'Ar', 'CO2'};
self.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'vP_vR', 0.5:0.01:2); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'SV');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
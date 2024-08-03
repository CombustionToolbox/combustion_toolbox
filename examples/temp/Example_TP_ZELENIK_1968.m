% -------------------------------------------------------------------------
% EXAMPLE: TP
%
% Compute equilibrium composition at defined temperature (e.g., 300 K) and
% pressure (e.g., 0.0101325 bar) for a stoichiometric CH2-O2 mixture, 
% and a set of 20 species considered
%   
% LS == {'C','CH','CH2','CH3','CH4','CO','CO2','C2','C2H2_acetylene',...
%        'C2H4','C2O','C3','H','HCO','HO2','H2','H2O','O','OH','O2'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update April 05 2022
% -------------------------------------------------------------------------

%% INITIALIZE
LS = {'C','CH','CH2','CH3','CH4','CO','CO2','C2','C2H2_acetylene',...
      'C2H4','C2O','C3','H','HCO','HO2','H2','H2O','O','OH','O2'};
self = App(LS);
%% INITIAL CONDITIONS
self = set_prop(self, 'pR', 0.01 * 1.01325, 'phi', 1);
self.PD.S_Fuel     = {'CH2'}; 
self.PD.S_Oxidizer = {'O2'}; 
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 300); 
%% TUNNUNG PARAMETERS
self.TN.tolN = 1e-14;
%% SOLVE PROBLEM
self = solve_problem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
post_results(self);